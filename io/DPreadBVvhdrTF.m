function [hdr] = DPread_brainvision_vhdrTF(filename); %input nSamplesPre added and later removed by DP...

%DP: it reads the header file IF and only IF they are exported in binary
%multiplexed IEEE float 32 data type

% READ_BRAINVISION_VHDR reads the known items from the BrainVision EEG
% header file and returns them in a structure
%
% Use as
%   hdr = read_brainvision_vhdr(filename)
%
% See also READ_BRAINVISION_EEG, READ_BRAINVISION_VMRK

% Copyright (C) 2003, Robert Oostenveld
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
% $Id: read_brainvision_vhdr.m 945 2010-04-21 17:41:20Z roboos $

hdr.DataFile         = read_asa(filename, 'DataFile=', '%s');
hdr.MarkerFile       = read_asa(filename, 'MarkerFile=', '%s');
hdr.DataFormat       = read_asa(filename, 'DataFormat=', '%s');
hdr.DataOrientation  = read_asa(filename, 'DataOrientation=', '%s');
hdr.BinaryFormat     = read_asa(filename, 'BinaryFormat=', '%s');
hdr.NumberOfChannels = read_asa(filename, 'NumberOfChannels=', '%d');
hdr.SamplingInterval = read_asa(filename, 'SamplingInterval=', '%f');   % microseconds

%----------------------------DP code---------------------------------------
hdr.DataType         = read_asa(filename, 'DataType=', '%s');
hdr.DataPoints       = read_asa(filename, 'DataPoints=', '%d');
hdr.SegmentationType = read_asa(filename, 'SegmentationType=', '%s');
hdr.SegmentDataPoints = read_asa(filename, 'SegmentDataPoints=', '%d');
hdr.Layers = read_asa(filename, 'Layers=', '%d');
hdr.LayerLowerLimit = read_asa(filename, 'LayerLowerLimit=', '%f');
hdr.LayerUpperLimit = read_asa(filename, 'LayerUpperLimit=', '%f');
hdr.LayerFunction = read_asa(filename, 'LayerFunction=', '%s');
hdr.Prop1 = read_asa(filename, 'Prop1=', '%s');
hdr.Prop2 = read_asa(filename, 'Prop2=', '%s');
hdr.Prop3 = read_asa(filename, 'Prop3=', '%s');
hdr.Prop4 = read_asa(filename, 'Prop4=', '%s');
hdr.Prop5 = read_asa(filename, 'Prop5=', '%s');
if strcmpi(hdr.DataType,'TIMEFREQUENCYDOMAIN_COMPLEX')
    hdr.ComplexDataPointFactor = 2;
else
    hdr.ComplexDataPointFactor = 1;
end
NumberOfcolumns = hdr.NumberOfChannels*hdr.Layers*hdr.ComplexDataPointFactor;
%----------------------------DP code---------------------------------------

if ~isempty(hdr.NumberOfChannels)
    for i=1:hdr.NumberOfChannels
        chan_str  = sprintf('Ch%d=', i);
        chan_info = read_asa(filename, chan_str, '%s');
        t = tokenize(chan_info, ',');
        
        if ~isempty(hdr.NumberOfChannels)
            for i=1:hdr.NumberOfChannels
                chan_str  = sprintf('Ch%d=', i);
                chan_info = read_asa(filename, chan_str, '%s');
                t = tokenize(chan_info, ',');
                
                if ~isempty(t)
                    
                    hdr.label{i} = t{1};
                    
                    if length(t)>1
                        hdr.reference{i}=t{2};
                    else
                        hdr.reference{i}='';
                    end
                    
                    hdr.resolution(i) = 1; %%1 is good for IEEE_FLOAT_32
                    if length(t)>2
                        tempRes = str2num(t{3});          % in microvolt
                        if isnumeric(tempRes)
                            hdr.resolution(i) = tempRes;
                        else
							hdr.resolution(i) = 1;
						end
                    end
                    
                    if length(t)>3
                        hdr.units{i}=t{4};
                    else
                        hdr.units{i}='';
                    end
                    
                else
                    hdr.label{i}='';
                    hdr.reference{i}='';
                    hdr.resolution(i) = 1;
                    hdr.units{i}='';
                end
                
                
                %----------------------------DP code---------------------------------------
                chan_info = DPread_asa(filename, chan_str, '%s',2);
                if ~isempty(chan_info)
                    t = tokenize(chan_info, ',');
                    hdr.coords.radius(i) = str2num(t{1});
                    hdr.coords.theta(i) = str2num(t{2});
                    hdr.coords.phi(i) = str2num(t{3});
                end
                %----------------------------DP code---------------------------------------
                
            end
        end
    end
end

% compute the sampling rate in Hz
hdr.Fs = 1e6/(hdr.SamplingInterval);

% the number of samples is unkown to start with
hdr.nSamples = Inf;

% determine the number of samples by looking at the binary file
if strcmpi(hdr.DataFormat, 'binary')
  % the data file is supposed to be located in the same directory as the header file
  % but that might be on another location than the present working directory
  [p, f, x] = fileparts(filename);
  datafile = fullfile(p, hdr.DataFile);
  info = dir(datafile);
  if isempty(info)
    error('cannot determine the location of the data file %s', hdr.DataFile);
  end
%   switch lower(hdr.BinaryFormat)
%     case 'int_16';
%       hdr.nSamples = info.bytes./(hdr.NumberOfChannels*2);
%     case 'int_32';
%      hdr.nSamples = info.bytes./(hdr.NumberOfChannels*4);
    if strcmpi(hdr.BinaryFormat,'ieee_float_32');
      hdr.nSamples = info.bytes./(NumberOfcolumns*4);
    else
        error('The data have to be exported from Vision Analyzer as ''IEEE_FLOAT_32'' ')
    end
%   end
% elseif strcmpi(hdr.DataFormat, 'ascii') && strcmpi(hdr.DataOrientation, 'vectorized')
%   % this is a very inefficient fileformat to read data from, it looks like this:
%   % Fp1   -2.129 -2.404 -18.646 -15.319 -4.081 -14.702 -23.590 -8.650 -3.957
%   % AF3   -24.023 -23.265 -30.677 -17.053 -24.889 -35.008 -21.444 -15.896 -12.050
%   % F7    -10.553 -10.288 -19.467 -15.278 -21.123 -25.066 -14.363 -10.774 -15.396
%   % F3    -28.696 -26.314 -35.005 -27.244 -31.401 -39.445 -30.411 -20.194 -16.488
%   % FC1   -35.627 -29.906 -38.013 -33.426 -40.532 -49.079 -38.047 -26.693 -22.852
%   % ...
%   fid = fopen(hdr.DataFile, 'rt');
%   tline = fgetl(fid);             % read the complete first line
%   fclose(fid);
%   t = tokenize(tline, ' ', true); % cut the line into pieces
%   hdr.nSamples = length(t) - 1;   % the first element is the channel label
else
    error('The data have to be exported from Vision Analyzer as ''binary'' ')
end

if isinf(hdr.nSamples)
  warning('cannot determine number of samples for this sub-fileformat');
end

%----------------------------DP code---------------------------------------
% the number of trials is unkown, assume continuous data
if ~isempty(hdr.SegmentDataPoints)
hdr.nTrials     = hdr.nSamples/hdr.SegmentDataPoints;
else
    hdr.nTrials=1;
end
%----------------------------DP code---------------------------------------

% ensure that the labels are in a column
hdr.label      = hdr.label(:);
hdr.reference  = hdr.reference(:);
hdr.resolution = hdr.resolution(:);
if all(isnan(hdr.resolution))
    hdr=rmfield(hdr,'resolution');
end



function [val] = DPread_asa(filename, elem, format, times, number, token)

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
    times=1;
end

if nargin<5
  if strcmp(format, '%s')
    number = 1;
  else
    number = Inf;
  end
end

if nargin<6
  token = '';
end


val = [];
elem = strtrim(lower(elem));

%----------------------------DP code modification--------------------------
STOP=0;
while STOP<times
  line = fgetl(fid);
  if ~isempty(line) && isequal(line, -1)
    % prematurely reached end of file
    fclose(fid);
    return
  end
  line = strtrim(line);
  lower_line = lower(line);
  if strmatch(elem, lower_line)
    STOP = STOP + 1;
  end
end
data = line((length(elem)+1):end);

while isempty(data)
  line = fgetl(fid);
  if isequal(line, -1)
    % prematurely reached end of file
    fclose(fid);
    return
  end
  data = strtrim(line);
end
%----------------------------DP code modification--------------------------

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

