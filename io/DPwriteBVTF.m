function DPwriteBVTF(dat, filename, hdr, absORpow, mrk) %, markerfile and absORpow DP addition, absORpow has to be 'abs' or 'pow'

% WRITE_BRAINVISION_EEG exports continuous EEG data to a BrainVision *.eeg
% and corresponding *.vhdr file. The samples in the exported file are
% multiplexed and stored in ieee-le float32 format.
%
% Use as
%   write_brainvision_eeg(filename, hdr, dat)
%
% See also READ_BRAINVISION_EEG, READ_BRAINVISION_VHDR, READ_BRAINVISION_VMRK

% Copyright (C) 2007, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: write_brainvision_eeg.m 4624 2011-10-29 10:10:49Z roboos $

%-------------------------------DP code------------------------------------
if strcmpi(absORpow,'abs')
    ABS='true';
    POW='false';
else
    POW='true';
    ABS='false';
    if isfield(hdr,'units')
        for i=1:hdr.NumberOfChannels; 
            hdr.units{i} = [hdr.units{i} '^2'];
        end
    end
end

% determine the filenames
[p, f, x] = fileparts(filename);
headerfile = fullfile(p, [f '.vhdr']);
datafile   = fullfile(p, [f ,x]);

if strcmpi(hdr.DataType,'TIMEFREQUENCYDOMAIN_COMPLEX')
    hdr.ComplexDataPointFactor = 2;
else
    hdr.ComplexDataPointFactor = 1;
end
NumberOfcolumns = hdr.NumberOfChannels*hdr.Layers*hdr.ComplexDataPointFactor;


disp('...restructuring data matrix...')
tic
nchan = size(dat,3);
nlrs  = size(dat,2);
nsmp = size(dat,1);
if length(size(dat))>3
    %DP: data comes as (time,frequency,channels,trials) and becomes
    %(trials,channels*frequency(*2),time)
    ntrl  = size(dat,4);
    dat=permute(dat,[2 3 1 4]); %(frequency,channels,time,trials)
    dat = reshape(dat,[nchan*nlrs,nsmp,ntrl]); %(frequency*channels,time,trials)
    if ~isreal(dat)
        datR=real(dat);
        datI=imag(dat);
        clear dat;
        dat(:,:,:,1)  = datR;
        clear datR;
        dat(:,:,:,2)  = datI; %(frequency*channels,time,trials,real/imag)
        clear datI;
        dat = permute(dat,[4,1,2,3]); %(real/imag,frequency*channels,time,trials)
        dat = reshape(dat,[NumberOfcolumns,nsmp,ntrl]); %(real/imag*frequency*channels,time,trials)
    end
    dat = permute(dat,[3,1,2]); %(trials,channels*frequency,time) or %(trials,real/imag*frequency*channels,time)
else
    %DP: data comes as (time,frequency,channels) and becomes
    %(channels*frequency(*2),time)
    dat=permute(dat,[2 3 1]); %(frequency,channels,time)
    dat = reshape(dat,[nchan*nlrs,nsmp]); %(frequency*channels,time)
    if ~isreal(dat)
        datR=real(dat);
        datI=imag(dat);
        clear dat;
        dat(:,:,1)  = datR; 
        clear datR;
        dat(:,:,2)  = datI; %(frequency*channels,time,real/imag)
        clear datI;
        dat = permute(dat,[3,1,2]); %(real/imag,frequency*channels,time)
        dat = reshape(dat,[NumberOfcolumns,nsmp]); %(real/imag*frequency*channels,time)
    end
end
toc

disp('...writing data file...')
tic
% open the data file and write the binary data
fid = fopen(datafile, 'wb', 'ieee-le');
if length(size(dat))>2
  %warning('writing segmented data as if it were continuous'); %DP commented this
  for i=1:ntrl
    fwrite(fid, squeeze(dat(i,:,:)), 'float32');
  end
else
  fwrite(fid, dat, 'float32');
end
fclose(fid);
toc
clear dat;



%-------------------------------DP code------------------------------------
disp('...writing marker file...')
if nargin>4
    
    if isstruct(mrk)
        markerfile = fullfile(p,[f,'.vmrk']);
        fid = fopen(markerfile, 'wb');
        fprintf(fid, 'Brain Vision Data Exchange Marker File, Version 2.0\r\n');
        fprintf(fid, '; Data created by FieldTrip modified by DP\r\n');
        fprintf(fid, '\r\n');
        fprintf(fid, '[Common Infos]\r\n');
        fprintf(fid, 'DataFile=%s\r\n',  [f ,x]);
        fprintf(fid, '\r\n');
        fprintf(fid, '[Marker Infos]\r\n');
%         ; Each entry: Mk<Marker number>=<Type>,<Description>,<Position in data points>,
%         ; <Size in data points>, <Channel number (0 = marker is related to all channels)>,
%         ; <Date (YYYYMMDDhhmmssuuuuuu)>, Visible
%         ; Fields are delimited by commas, some fields may be omitted (empty).
%         ; Commas in type or description text are coded as "\1".
        Nmarkers = length(mrk.Position);
        for i=1:Nmarkers;
            if strcmpi(mrk.type{i},'New Segment')
                fprintf(fid, 'Mk%d=%s,%s,%d,%d,%d,%s\r\n', i, mrk.type{i}, mrk.label{i},mrk.Position(i),mrk.Points(i),mrk.Channel(i),mrk.Date{i});
            else
                fprintf(fid, 'Mk%d=%s,%s,%d,%d,%d\r\n', i, mrk.type{i}, mrk.label{i},mrk.Position(i),mrk.Points(i),mrk.Channel(i));
            end
        end
        fprintf(fid, '\r\n');
        fprintf(fid, '[Marker User Infos]\r\n');
        fclose(fid);
    elseif ischar(mrk)
        if exist(mrk,'file')
            markerfile = fullfile(p,[f,'.vmrk']);
            markerfileOLD=mrk;
            copyfile(markerfileOLD,markerfile,'f');
            MarkerDataFile = read_asa(markerfile, 'DataFile=', '%s');
            findreplace(markerfile,['DataFile=',MarkerDataFile],['DataFile=',[f '.eeg']])
        else
            markerfile = '';
        end
        
    else
        markerfile = '';
    end
    
else
    markerfile = '';
end
%-------------------------------DP code------------------------------------




disp('...writing header file...')

if hdr.NumberOfChannels~=nchan
  error('number of channels in header does not match with the data');
end

% this is the only supported data format
hdr.DataFormat      = 'BINARY';
hdr.DataOrientation = 'MULTIPLEXED';
hdr.BinaryFormat    = 'IEEE_FLOAT_32';
hdr.resolution      = ones(hdr.NumberOfChannels,1);  % no additional calibration needed, since float32

% open the header file and write the ascii header information
fid = fopen(headerfile, 'wb');
fprintf(fid, 'Brain Vision Data Exchange Header File Version 2.0\r\n');
fprintf(fid, '; Data created by FieldTrip modified by DP\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '[Common Infos]\r\n');
fprintf(fid, 'DataFile=%s\r\n',          [f ,x]);
if ~isempty(markerfile) 
  fprintf(fid, 'MarkerFile=%s\r\n',      [f , '.vmrk']);
end
fprintf(fid, 'DataFormat=%s\r\n',        'BINARY'); %DP changed it from hdr.DataFormat
fprintf(fid, 'DataOrientation=%s\r\n',   'MULTIPLEXED');  %DP changed it from hdr.DataOrientation
fprintf(fid, 'DataType=%s\r\n',   hdr.DataType); 
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

%--------------------these are special for TF data--------=====------------
fprintf(fid, 'Layers=%d\r\n',  hdr.Layers);
fprintf(fid, 'LayerLowerLimit=%f\r\n',  hdr.LayerLowerLimit);
fprintf(fid, 'LayerUpperLimit=%f\r\n',  hdr.LayerUpperLimit);
fprintf(fid, 'LayerFunction=%s\r\n',   hdr.LayerFunction);
fprintf(fid, '\r\n');
fprintf(fid, '[User Infos]\r\n');
if  ~isempty(hdr.Prop1)
    fprintf(fid, 'Prop1=%s\r\n',   hdr.Prop1);
else
    fprintf(fid, 'Prop1=%s\r\n',   ['int,BrainVision.WTType,2']);%????TEST?????
end
if  ~isempty(hdr.Prop2)
    fprintf(fid, 'Prop2=%s\r\n',   hdr.Prop2);
else
    fprintf(fid, 'Prop2=%s\r\n',   ['bool,BrainVision.WAbsVal,',ABS]);%????TEST?????
end
if ~isempty(hdr.Prop3)
    fprintf(fid, 'Prop3=%s\r\n',   hdr.Prop3);
else
    fprintf(fid, 'Prop3=%s\r\n',   ['bool,BrainVision.WPowVal,',POW]);%????TEST?????
end
if  ~isempty(hdr.Prop4)
    fprintf(fid, 'Prop4=%s\r\n',   hdr.Prop4);
else
    fprintf(fid, 'Prop4=%s\r\n',   ['string,BrainVision.CWTWavelet,Morlet Complex']);%????TEST?????
end
if isfield(hdr,'Prop5') && ~isempty(hdr.Prop5)
    fprintf(fid, 'Prop5=%s\r\n',   hdr.Prop5);
else
    fprintf(fid, 'Prop5=%s\r\n',   ['single,BrainVision.MorletFactor,3']);%????TEST?????
end
%--------------------these are special for TF data--------=====------------

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







function [val] = read_asa(filename, elem, format, number, token)

% READ_ASA reads a specified element from an ASA file
%
% val = read_asa(filename, element, type, number)
%
% where the element is a string such as
%   NumberSlices
%   NumberPositions
%   Rows
%   Columns
%   etc.
%
% and format specifies the datatype according to
%   %d  (integer value)
%   %f  (floating point value)
%   %s  (string)
%
% number is optional to specify how many lines of data should be read
% The default is 1 for strings and Inf for numbers.
%
% token is optional to specifiy a character that separates the values from
% anything not wanted.

% Copyright (C) 2002, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: read_asa.m 945 2010-04-21 17:41:20Z roboos $

fid = fopen(filename, 'rt');
if fid==-1
  error(sprintf('could not open file %s', filename));
end

if nargin<4
  if strcmp(format, '%s')
    number = 1;
  else
    number = Inf;
  end
end

if nargin<5
  token = '';
end


val = [];
elem = strtrim(lower(elem));

while (1)
  line = fgetl(fid);
  if ~isempty(line) && isequal(line, -1)
    % prematurely reached end of file
    fclose(fid);
    return
  end
  line = strtrim(line);
  lower_line = lower(line);
  if strmatch(elem, lower_line)
    data = line((length(elem)+1):end);
    break
  end
end

while isempty(data)
  line = fgetl(fid);
  if isequal(line, -1)
    % prematurely reached end of file
    fclose(fid);
    return
  end
  data = strtrim(line);
end

if strcmp(format, '%s')
  if number==1
    % interpret the data as a single string, create char-array
    val = detoken(strtrim(data), token);
    fclose(fid);
    return
  end
  % interpret the data as a single string, create cell-array
  val{1} = detoken(strtrim(data), token);
  count = 1;
  % read the remaining strings
  while count<number
    line = fgetl(fid);
    if ~isempty(line) && isequal(line, -1)
      fclose(fid);
      return
    end
    tmp = sscanf(line, format);
    if isempty(tmp)
      fclose(fid);
      return
    else
      count = count + 1;
      val{count} = detoken(strtrim(line), token);
    end
  end

else
  % interpret the data as numeric, create numeric array
  count = 1;
  data = sscanf(detoken(data, token), format)';
  if isempty(data),
    fclose(fid);
    return
  else
    val(count,:) = data;
  end
  % read remaining numeric data
  while count<number
    line = fgetl(fid);
    if ~isempty(line) && isequal(line, -1)
      fclose(fid);
      return
    end
    data = sscanf(detoken(line, token), format)';
    if isempty(data)
      fclose(fid);
      return
    else
      count = count+1;
      val(count,:) = data;
    end
  end
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = detoken(in, token)
if isempty(token)
  out = in;
  return;
end

[tok rem] = strtok(in, token);
if isempty(rem)
  out = in;
  return;
else
  out = strtok(rem, token);
  return
end



function [tok] = tokenize(str, sep, rep)

% TOKENIZE cuts a string into pieces, returning the pieces in a cell array
%
% Use as
%   t = tokenize(str)
%   t = tokenize(str, sep)
%   t = tokenize(str, sep, rep)
% where
%   str = the string that you want to cut into pieces
%   sep = the separator at which to cut (default is whitespace)
%   rep = whether to treat repeating seperator characters as one (default is false)
%
% With the optional boolean flag "rep" you can specify whether repeated
% seperator characters should be squeezed together (e.g. multiple
% spaces between two words). The default is rep=1, i.e. repeated
% seperators are treated as one.
%
% See also STRTOK, TEXTSCAN

% Copyright (C) 2003-2010, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: tokenize.m 4624 2011-10-29 10:10:49Z roboos $

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

if nargin<2
  sep = [9:13 32]; % White space characters
end

if nargin<3
  rep = false;
end

current_argin = {str, sep, rep};
if isequal(current_argin, previous_argin)
  % don't do the processing again, but return the previous values from cache
  tok = previous_argout;
  return
end

tok = {};
f = find(ismember(str, sep));
f = [0, f, length(str)+1];
for i=1:(length(f)-1)
  tok{i} = str((f(i)+1):(f(i+1)-1));
end

if rep
  % remove empty cells, which occur if the separator is repeated (e.g. multiple spaces)
  tok(cellfun('isempty', tok))=[];
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout  = tok;
previous_argin  = current_argin;
previous_argout = current_argout;


function findreplace(file,otext,ntext,varargin)

%FINDREPLACE finds and replaces strings in a text file
%
% SYNTAX:
%
% findreplace(file,otext,ntext)
% findreplace(file,otext,ntext,match)
%
% findreplace  : This function finds and replaces strings in a text file
%
%       file:           text file name (with or without path)
%       otext:          text to be replaced (old text)
%       ntext:          replacing text (new text)
%       match:          either (1) for match case or (0) to ignore case
%                       default value is (1)
%
% Example:
%   findreplace('sample.txt','Moller','Moler');
%   findreplace('sample.txt','jake','Jack',0);
%   findreplace('sample.txt','continue it is','continue its',0);
%
%   Copyright 2005 Fahad Al Mahmood
%   Version: 1.0    $Date: 24-Dec-2005

% Obtaining the file full path
[fpath,fname,fext] = fileparts(file);
if isempty(fpath)
    out_path = pwd;
elseif fpath(1)=='.'
    out_path = [pwd filesep fpath];
else
    out_path = fpath;
end

% Reading the file contents
k=1;
all=0;
opt=[];
first_time=1;
change_counter=0;
fid = fopen([out_path filesep fname fext],'r');
while 1
    line{k} = fgetl(fid);
    if ~ischar(line{k})
        break;
    end
    k=k+1;
end
fclose(fid);
old_lines = line;

%Number of lines
nlines = length(line)-1;

for i=1:nlines
    if nargin==3, match=1;
    else match=varargin{1}; end

    if match==1, loc = regexp(line{i},otext);
    elseif match==0, loc = regexpi(line{i},otext);
    end

    if ~isempty(loc)
        nloc = 1;
        for j=1:length(loc)
            if all==0    
%-------------------------------DP code------------------------------------
                % Displaying keyboard instructions
%                 if first_time
%                     disp(' ');
%                     disp('(y) change (n) skip (a) change all (s) stop');
%                     disp(' ');
%                     first_time=0;
%                 end
%                 disp(line{i});
%                 opt = input(underline(loc(j),length(otext)),'s');

%                if opt=='a'
                    line{i} = regexprep(line{i},otext,ntext, 'preservecase',nloc);
                    change_counter = change_counter + 1;
                    all=1;
                    if length(loc)>j
                        loc(j:end) = loc(j:end) + (length(ntext)-length(otext));
                    end
%                 elseif opt=='y'
%                     line{i} = regexprep(line{i},otext,ntext, 'preservecase',nloc);
%                     change_counter = change_counter + 1;
%                     if length(loc)>j
%                         loc(j:end) = loc(j:end) + (length(ntext)-length(otext));
%                     end
%                 elseif opt=='s';
%                     break;
%                 else
%                     nloc = nloc + 1;
%                 end
%-------------------------------DP code------------------------------------
            else
                line{i} = regexprep(line{i},otext,ntext, 'preservecase',nloc);
                change_counter = change_counter + 1;
                if length(loc)>j
                    loc(j:end) = loc(j:end) + (length(ntext)-length(otext));
                end
            end
        end
    end
    if opt=='s';
        break
    end
end

line = line(1:end-1);

%disp(' ');%DP commented this
if change_counter~=0
    % Writing to file
    fid2 = fopen([out_path filesep fname fext],'w');

    for i=1:nlines
        fprintf(fid2,[line{i} '\n']);
    end
    fclose(fid2);
    %disp([num2str(change_counter) ' Changes Made & Saved Successfully']); %DP commented this
    %disp(' '); %DP commented this
else
    disp('No Match Found / No Change Applied');
    disp(' ');
end

function uline = underline(loc,length)
s=' ';
l='-';
uline=[];
if loc==1
    for i=1:length
        uline=[uline l];
    end
else
    for i=1:loc-1
        uline = [uline s];
    end
    for i=1:length
        uline=[uline l];
    end
end
uline = [uline s];







