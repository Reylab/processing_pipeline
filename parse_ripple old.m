function parse_ripple(filenames,remove_chs,macro,max_memo_GB, overwrite, channels)

% This code requires the neuroshape library in the path.
% max_memo_GB is an idea of the number of GB allocated for the data to be
% stored in RAM, so it is used to compute the number of segments in which
% the data should be split for processing

if ~exist('filenames','var')
   aux=dir('*.ns5');
   filenames= {aux.name};
end

expr2remove = '-\d+$';
%%
if ~isempty(mfilename)
    root_rc = [fileparts(mfilename('fullpath')) filesep '..' filesep '..'];
    if exist([root_rc filesep 'reylab_custompath.m'],'file')
        addpath(root_rc);
        custom_path = reylab_custompath('neuroshare');
    end
end

%%


if ~exist('overwrite','var') || isempty(overwrite)
    overwrite = false;
end

if ~exist('macro','var')
    macro = [];
end

if ~exist('remove_chs','var')
    remove_chs = [];
end
with_memory=true;

try
	memory;
catch
	with_memory=false;
end
if with_memory
	[~, systemview] = memory;
	memo_avaible = floor(systemview.PhysicalMemory.Available*0.80);
	if exist('max_memo_GB','var') && ~isempty(max_memo_GB)
        max_memo = max_memo_GB*(1024)^3;
		if max_memo > memo_avaible
			error('max_memo_GB > 80% of Physical Memory Available')
		end
	else
		max_memo = memo_avaible;
	end
else
    if ~exist('max_memo_GB','var') || isempty(max_memo_GB)
        warning('max_memo_GB must be set. Setting max memo as 8GB')
        max_memo_GB = 8;
	end
	max_memo = max_memo_GB*(1024)^3;
end

tcum=0;

if ischar(filenames)
    filenames = {filenames};
end

formatvector=@(v,f) sprintf(['[' repmat(['%' f ', '], 1, numel(v)-1), ['%' f ']\n']  ],v);

metadata_file = fullfile(pwd, 'NSx.mat');
if exist(metadata_file,'file')
    metadata = load(metadata_file);
else
    metadata = [];
end

for fi = 1:length(filenames)
    filename= filenames{fi};
    new_files(fi).name = filename;
    if length(filename)<3 || (~strcmpi(filename(2:3),':\') && ...
                     ~strcmpi(filename(1),'/') && ...
                     ~strcmpi(filename(2),'/') && ...
                     ~strcmpi(filename(1:2), '\\')&& ~strcmpi(filename(2:3),':/'))
        filename= [pwd filesep filename];
    end
    
    tic
    [ns_status, hFile] = ns_OpenFile(filename, 'single');
    if strcmp(ns_status,'ns_FILEERROR')
        error('Unable to open file: %s',filename)
    end
    
    switch hFile.FileInfo.Type(3)
        case '1'
            sr = 500;
        case '2'
            sr = 1000;
        case '3'
            sr = 2000;
        case '4'
            sr = 10000;            
        case {'5' ,'6'}
            sr = 30000;
        otherwise
            error('ERROR: %s file type not supported',hFile.FileInfo.Type)
    end
    
    nchan = size(hFile.Entity,2);   % number of channels
    samples_per_channel = ceil(max_memo/(nchan*2));
    if fi == 1
        %get info about channels in nsx file
        chs_info = struct();
        chs_info.unit = cellfun(@(x) deblank(x),{hFile.Entity.Units},'UniformOutput',false);
        chs_info.label = cellfun(@(x) deblank(x),{hFile.Entity.Label},'UniformOutput',false)';
        chs_info.conversion = cell2mat({hFile.Entity.Scale});
        chs_info.id = cell2mat({hFile.Entity.ElectrodeID});
        chs_info.pak_list = 0*chs_info.id;
        
        chs_info.macro = chs_info.label;
        micros = cellfun(@(x) strcmp(x,'uV'),chs_info.unit);
        if ~isempty(macro)
            chs_info.macro(micros) = arrayfun(@(x) macro{ceil(x/9)},find(micros),'UniformOutput',false);
        end
        
        outfile_handles = cell(1,nchan); %some will be empty
        [~,~,fext] = fileparts(filename);
        fext = lower(fext(2:end));
        nsx_ext = fext(end);
        ch_ext = ['.NC' nsx_ext];
        if ~exist('channels','var') || isempty(channels)
            channels = hFile.FileInfo.ElectrodeList;
        end
        if ~isempty(remove_chs)
            channels = setdiff(channels, remove_chs);
        end
        parsed_chs = [];
        new_channel_id = [];
        for i = 1:nchan
            c = hFile.FileInfo.ElectrodeList(i);
            if ismember(c,channels)
                ccname = c;
                if ~isempty(metadata) %NSx file in current folder
                    repetead = arrayfun(@(x) (x.chan_ID==ccname) && (x.sr==sr) ,metadata.NSx);
                    if ~isempty(repetead) && sum(repetead)>0 %found channel
                        pos = find(repetead);
                        if overwrite
                            f2delete = [metadata.NSx(pos).output_name  metadata.NSx(pos).ext];
                            fprintf('Overwritting channel %d. Deleting file %s\n', metadata.NSx(pos).chan_ID, f2delete)
                            delete(f2delete)
                        else
                            fprintf('Skipping channel %d, already parsed.\n', metadata.NSx(pos).chan_ID)
                            continue %If output_name wasn't set, the existing parsed channels won't be overwritten.
                        end
                    end
                end
                parsed_chs(end+1) = c;
                new_channel_id(end+1) = ccname;

                ix = chs_info.id==c;
                %output_name = chs_info.label{ix};
                output_name = chs_info.macro{ix};
                outn_i = regexp(output_name,expr2remove);
                if ~isempty(outn_i) && outn_i>1
                    output_name = output_name(1:outn_i-1);
                end
                chs_info.output_name{ix} = [output_name '_' num2str(ccname)];
                outfile_handles{i} = fopen([chs_info.output_name{ix} ch_ext],'w');
            end
        end

        new_files(fi).first_sample = 1;
    else
        new_files(fi).first_sample = new_files(fi-1).lts + new_files(fi-1).first_sample;
    end
    if isempty(parsed_chs)
        disp('Without channels to parse.')
        return
    end
   
    %total lenght adding the zeros from Timestamp
    lts = hFile.TimeSpan;
    new_files(fi).lts = lts;

    N = lts;   % total data points
    num_segments = ceil(N/samples_per_channel);
    fprintf('Data will be processed in %d segments of %d samples each.\n',num_segments,min(samples_per_channel,N))
    
    min_valid_val = zeros([nchan,1]);
    for i = 1:nchan
        [~, nsAnalogInfo] = ns_GetAnalogInfo(hFile, i);
        min_valid_val(i) = nsAnalogInfo.MinVal/nsAnalogInfo.Resolution;
    end
    
    for j=1:num_segments
        ini = (j-1)*samples_per_channel+1;
        fin = min(j*samples_per_channel,N);
        tcum = tcum + toc;  % this is because openNSx has a tic at the beginning
        [~, Data] = ns_GetAnalogDataBlock(hFile, [1:nchan], ini, fin-ini+1, 'unscale');
        for i = 1:nchan
            if ~isempty(outfile_handles{i}) %channels with empty outfile_handles{i} are not selected
                pak_lost = Data(1:end,i)<min_valid_val(i);
                Data(pak_lost,i)=0;
                chs_info.pak_list(i) = chs_info.pak_list(i)+sum(pak_lost);
                fwrite(outfile_handles{i},Data(1:end,i),'int16');
            end
        end
        fprintf('Segment %d out of %d processed. ',j,num_segments)
    end
    
    tcum = tcum + toc;
    fprintf('Total time spent in parsing the data was %s secs.\n',num2str(tcum, '%0.1f')); 
end
fclose('all');



if isempty(metadata)
    files = [];
    NSx = [];
else
    NSx = metadata.NSx;
    files = metadata.files;
end
lts = sum([new_files(:).lts]);
fprintf('%d data points written per channel\n',lts)


for ci = 1:length(new_channel_id)
    ch = new_channel_id(ci);
    elec_id = parsed_chs(ci);
    repetead = arrayfun(@(x) (x.chan_ID==ch) && (x.sr==sr) ,NSx);
    if isempty(repetead) || sum(repetead)==0
        pos = length(NSx)+1;
    else
        pos = find(repetead);
    end
        ix = chs_info.id==elec_id;
        NSx(pos).chan_ID = ch;
        NSx(pos).conversion = chs_info.conversion(ix);
        NSx(pos).label = chs_info.label{ix};
        NSx(pos).macro = get_macro(chs_info.macro{ix});
        NSx(pos).unit = chs_info.unit{ix};
        NSx(pos).electrode_ID = elec_id;
        NSx(pos).which_system = 'RIPPLE';
        NSx(pos).ext = ch_ext;
        NSx(pos).lts = lts;
        NSx(pos).filename = filenames;
        NSx(pos).sr = sr;
        NSx(pos).output_name = chs_info.output_name{ix};
end

for i = 1:length(new_files)
    repetead = arrayfun(@(x) strcmp(x.name,new_files(i).name),files);
    if isempty(repetead) || sum(repetead)==0
        pos = length(files)+1;
    else
        pos = find(repetead);
    end
    files(pos).name = new_files(i).name;
    files(pos).first_sample = new_files(i).first_sample;
    files(pos).lts = new_files(i).lts;
end
freq_priority=[30000, 2000, 10000, 1000, 500];
save(metadata_file, 'NSx','files','freq_priority')
custom_path.rm()
end

function macro=get_macro(label)
    macro = regexp(label,'^\w*(?=(\d* raw$))','match');
    if ~isempty(macro)
        macro = macro{1};
    else
        macro = label;
    end
end
