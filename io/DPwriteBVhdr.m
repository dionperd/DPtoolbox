function DPwriteBVhdr(filename, hdr, markerfile)


if nargin<3
    markerfile='';
end

% determine the filenames
[p, f, x] = fileparts(filename);
headerfile = fullfile(p, [f '.vhdr']);

disp('...writing header file...')

% this is the only supported data format
hdr.DataFormat      = 'BINARY';
hdr.DataOrientation = 'MULTIPLEXED';
hdr.BinaryFormat    = 'IEEE_FLOAT_32';
hdr.resolution      = ones(hdr.NumberOfChannels,1);  % no additional calibration needed, since float32
% open the header file and write the ascii header information
fid = fopen(headerfile, 'wb');
fprintf(fid, 'Brain Vision Data Exchange Header File, Version 2.0\r\n');
fprintf(fid, '; Data created by FieldTrip modified by DP\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '[Common Infos]\r\n');
fprintf(fid, 'DataFile=%s\r\n',          [f , x]);
if ~isempty(markerfile) 
  fprintf(fid, 'MarkerFile=%s\r\n',      [f , '.vmrk']);
end
fprintf(fid, 'DataFormat=%s\r\n',        'BINARY'); %DP changed it from hdr.DataFormat
fprintf(fid, 'DataOrientation=%s\r\n',   'MULTIPLEXED');  %DP changed it from hdr.DataOrientation
fprintf(fid, 'DataType=%s\r\n',   'TIMEDOMAIN');
fprintf(fid, 'NumberOfChannels=%d\r\n',  hdr.NumberOfChannels);
fprintf(fid, 'DataPoints=%d\r\n',  hdr.DataPoints);
% Sampling interval in microseconds
%-------------------------------DP code------------------------------------
fprintf(fid, 'SamplingInterval=%d\r\n',  hdr.SamplingInterval);
if isfield(hdr,'SegmentationType')
    fprintf(fid, 'SegmentationType=%s\r\n',   hdr.SegmentationType);
end
if isfield(hdr,'SegmentDataPoints')
    fprintf(fid, 'SegmentDataPoints=%d\r\n',   hdr.SegmentDataPoints);
end
%-------------------------------DP code------------------------------------
fprintf(fid, '\r\n');
fprintf(fid, '[Binary Infos]\r\n');
fprintf(fid, 'BinaryFormat=%s\r\n',      'IEEE_FLOAT_32'); %DP changed it from hdr.BinaryFormat
fprintf(fid, '\r\n');
fprintf(fid, '[Channel Infos]\r\n');
% Each entry: Ch<Channel number>=<Name>,<Reference channel name>,<Resolution in microvolts>,<Future extensions>...
% Fields are delimited by commas, some fields might be omitted (empty).
% Commas in channel names should be coded as "\1", but are not supported
% here
if ~isfield(hdr,'resolution')
    for i=1:hdr.NumberOfChannels; 
        hdr.resolution{i} = '';
    end
else
    for i=1:hdr.NumberOfChannels;
        if isnan(hdr.resolution(i))
            temp{i} = '';
        else
            temp{i} = num2str(hdr.resolution(i));
        end
    end
    hdr.resolution=temp;
    clear temp;
end
        
if ~isfield(hdr,'units')
    for i=1:hdr.NumberOfChannels; 
        hdr.units{i} = '';
    end
end

if ~isfield(hdr,'reference')
    for i=1:hdr.NumberOfChannels; 
        hdr.reference{i} = '';
    end
end

if ~isfield(hdr,'label')
    for i=1:hdr.NumberOfChannels; 
        hdr.label{i} = '';
    end
end

for i=1:hdr.NumberOfChannels;
    fprintf(fid, 'Ch%d=%s,%s,%s,%s\r\n', i, hdr.label{i},hdr.reference{i}, hdr.resolution{i},hdr.units{i});
end


if isfield(hdr,'coords')
    if isfield(hdr.coords,'radius') && isfield(hdr.coords,'theta') && isfield(hdr.coords,'phi')
        if (numel(hdr.coords.radius) == hdr.NumberOfChannels) && (numel(hdr.coords.theta) == hdr.NumberOfChannels) && (numel(hdr.coords.phi) == hdr.NumberOfChannels)
            fprintf(fid, '\r\n');
            fprintf(fid, '[Coordinates]\r\n');
            for i=1:hdr.NumberOfChannels %DP changed it from hdr.NumberOfChannels
                %; Each entry: Ch<Channel number>=<Radius>,<Theta>,<Phi>
                fprintf(fid, 'Ch%d=%d,%d,%g\r\n', i, hdr.coords.radius(i), hdr.coords.theta(i), hdr.coords.phi(i));
            end
        end
    end
end
fclose(fid);