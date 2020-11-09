function parsed_chs = parse_data_NSx(filenames,max_memo_GB,output_name,channels)
% Code only valid for recordings without pauses or data loss (Nsx.Data field can't be a cell)
% This code requires the file openNSx.m, from the NPMK, in the path. You can download the download NPMK from https://github.com/BlackrockMicrosystems/NPMK/releases
% max_memo_GB is an idea of the number of GB allocated for the data to be
% stored in RAM, so it is used to compute the number of segments in which
% the data should be split for processing
if ~exist('output_name','var') || isempty(output_name)
    output_name = 'NSX';
end


with_memory=true;
try
	memory;
catch
	with_memory=false;
end
if with_memory
	[userview,systemview] = memory;
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
	max_memo = max_memo_GB*(1024)^3;
end

tcum=0;
if ischar(filenames)
    filenames = {filenames};
end

for fi = 1:length(filenames)
    filename= filenames{fi};
    files(fi).name = filename;
    if length(filename)<3 || (~strcmpi(filename(2:3),':\') && ...
                     ~strcmpi(filename(1),'/') && ...
                     ~strcmpi(filename(2),'/') && ...
                     ~strcmpi(filename(1:2), '\\')&& ~strcmpi(filename(2:3),':/'))

        filename= [pwd filesep filename];
    end
    
    NSx = openNSx(filename, 'report','noread');
    nchan = NSx.MetaTags.ChannelCount;   % number of channels
    sr = NSx.MetaTags.SamplingFreq;   % sampling rate
    lts = sum(NSx.MetaTags.DataPoints);   % total data points 
    files(fi).lts = lts;
    if fi == 1
        outfile_handles = cell(1,nchan); %some will be empty
        [~,~,fext] = fileparts(filename);
        fext = lower(fext(2:end));
        nsx_ext = fext(end);
        ch_ext = ['.NC' nsx_ext];
        chs_info = table();
        chs_info.unit = cellfun(@(x) max_analog2unit(x),{NSx.ElectrodesInfo.MaxAnalogValue},'UniformOutput',false)';
        chs_info.label = cellfun(@(x) deblank(x),{NSx.ElectrodesInfo.Label},'UniformOutput',false)';
        chs_info.conversion = (double(cell2mat({NSx.ElectrodesInfo.MaxAnalogValue}))./double(cell2mat({NSx.ElectrodesInfo.MaxDigiValue})))';
        chs_info.Properties.RowNames = cellfun(@(x) num2str(x),{NSx.ElectrodesInfo.ElectrodeID},'UniformOutput',false)';
        if ~exist('channels','var') || isempty(channels)
            channels = NSx.MetaTags.ChannelID;
        end
        parsed_chs = [];
        for i = 1:nchan
            c = NSx.MetaTags.ChannelID(i);
            if ismember(c,channels)
                parsed_chs(end+1) = c;
                outfile_handles{i} = fopen([output_name '_' num2str(c) ch_ext],'w');
            end
        end
        files(fi).first_sample = 1;
    else
        files(fi).first_sample = files(fi-1).lts + files(fi-1).first_sample;
    end
    

    DataPoints = NSx.MetaTags.DataPoints;

    DataPoints = [0 DataPoints]; %added for use as starting index
    samples_per_channel = ceil(max_memo/(nchan*length(NSx.MetaTags.DataPoints))/2);
    for part = 2:length(DataPoints)
        N = DataPoints(part);   % total data points 
        num_segments = ceil(N/samples_per_channel);
        fprintf('TimeStamps generated. Data will be processed in %d segments of %d samples each.\n',num_segments,min(samples_per_channel,N))
        for j=1:num_segments
            ini = (j-1)*samples_per_channel+1+DataPoints(part-1);
            fin = min(j*samples_per_channel,N)+DataPoints(part-1);
            tcum = tcum + toc;  % this is because openNSx has a tic at the beginning
            NSx = openNSx('read',filename,['t:' num2str(ini) ':' num2str(fin)]);
            fix_ini = 0;
            if size(NSx.Data,2) > (fin-ini +1)
                fix_ini = size(NSx.Data,2) - (fin-ini +1);
                warning(sprintf('segment %d: %d extra samples will be removed.',j,fix_ini))
            end
            for i = 1:nchan
                if ~isempty(outfile_handles{i})
                    %BE AWARE THAT THE RAW DATA IN THE NS5 AND NC5 IS SCALED UP BY A FACTOR OF 4
                    fwrite(outfile_handles{i},NSx.Data(i,(1+fix_ini):end),'int16');
                end
            end
            fprintf('Segment %d out of %d processed. Data Point Read = %d \n',j,num_segments,size(NSx.Data,2));
        end
    end
    
    tcum = tcum + toc;
    fprintf('Total time spent in parsing the data was %s secs.\n',num2str(tcum, '%0.1f')); 
end
fclose('all');
metadata_file = fullfile(pwd, [output_name '_TimeStamps.mat']);
metadata = struct;
if exist(metadata_file,'file')
    metadata = load(metadata_file);
end
lts = sum([files(:).lts]);
fprintf('%d data points written per channel\n',lts)
if strcmp(nsx_ext,'5') 
    if isfield(metadata,'parsed_chs')
        parsed_chs = union(parsed_chs,metadata.parsed_chs);
    end
    metadata.parsed_chs = parsed_chs;
    metadata.lts = lts;
    metadata.nchan = nchan;
    metadata.sr = sr;
    metadata.ch_ext = ch_ext;
    metadata.chs_info = chs_info;
    metadata.files = files;
    save(metadata_file, '-struct', 'metadata')
else
    if isfield(metadata,fext) && isfield(metadata.(fext),'sr')
        parsed_chs = union(parsed_chs,metadata.(fext).parsed_chs);
    end
    metadata.(fext).parsed_chs = parsed_chs;
    metadata.(fext).lts = lts;
    metadata.(fext).nchan = nchan;
    metadata.(fext).sr = sr;
    metadata.(fext).ch_ext = ch_ext;
    metadata.(fext).chs_info = chs_info;
    metadata.(fext).files = files;
    save(metadata_file, '-struct', 'metadata')
end
end
function unit = max_analog2unit(x)
    switch x
        case 5000
            unit='mV';
        case 8191
            unit ='uV';
    end
end
