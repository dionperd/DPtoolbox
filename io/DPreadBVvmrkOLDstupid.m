function allinfo = DPreadBVvmrkOLDstupid(filename)
% READ_BRAINVISION_VMRK reads the markers and latencies
% it returns the stimulus/response code and latency in ms.
%
% Use as
%   [stim, resp, segment, timezero] = read_brainvision_vmrk(filename)
% 
% This function needs to read the header from a separate file and
% assumes that it is located at the same location.
%
% See also READ_BRAINVISION_VHDR, READ_BRAINVISION_EEG

% original M. Schulte 31.07.2003
% modifications R. Oostenveld 14.08.2003
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
% $Id: read_brainvision_vmrk.m 945 2010-04-21 17:41:20Z roboos $

% stimulus={};
% response={};
% segment=[];
% timezero=[];
allinfo=struct();

fid=fopen(filename,'rt');
if fid==-1,
    error('cannot open marker file')
end

iM=0;
line=1;
while line~=-1,
    line=fgetl(fid);
    % pause
    if ~isempty(line),
        if ~isempty(findstr(line,'Mk')),
            if ~isempty(findstr(line,'Stimulus'))
                iM=iM+1;
                allinfo.type{iM} = 'Stimulus';
                [token,rem] = strtok(line,',');
                [type,rem] = strtok(rem,',');
                allinfo.label{iM}= type;
                [token,rem] = strtok(rem,',');
                allinfo.Position(iM) = sscanf(token,'%i');
                %stimulus=[stimulus; {type allinfo.Position(iM)}];
                [token,rem] = strtok(rem,',');
                allinfo.Points(iM) = str2num(token);
                [token,rem] = strtok(rem,',');
                allinfo.Channel(iM) = str2num(token);
                allinfo.Date{iM} = '';
                
            elseif ~isempty(findstr(line,'Response'))
                iM=iM+1;
                allinfo.type{iM} = 'Response';
                [token,rem] = strtok(line,',');
                [type,rem] = strtok(rem,',');
                allinfo.label{iM}= type;
                [token,rem] = strtok(rem,',');
                allinfo.Position(iM) = sscanf(token,'%i');
                %response=[response; {type allinfo.Position(iM)}];
                [token,rem] = strtok(rem,',');
                allinfo.Points(iM) = str2num(token);
                [token,rem] = strtok(rem,',');
                allinfo.Channel(iM) = str2num(token);
                allinfo.Date{iM} = '';
                
            elseif ~isempty(findstr(line,'New Segment'))
                iM=iM+1;
                allinfo.type{iM} = 'New Segment';
                [token,rem] = strtok(line,',');
                allinfo.label{iM} =nan;
                [token,rem] = strtok(rem,',');
                allinfo.Position(iM) = sscanf(token,'%i');
                %segment=[segment; allinfo.Position(iM)];
                [token,rem] = strtok(rem,',');
                allinfo.Points(iM) = str2num(token);
                [token,rem] = strtok(rem,',');
                allinfo.Channel(iM) = str2num(token);
                [token,rem] = strtok(rem,',');
                allinfo.Date{iM} = token;
                
            elseif ~isempty(findstr(line,'Time 0'))
                iM=iM+1;
                allinfo.type{iM} = 'Time 0';
                [token,rem] = strtok(line,',');
                allinfo.label{iM} =nan;
                [token,rem] = strtok(rem,',');
                allinfo.Position(iM) = sscanf(token,'%i');
                %timezero=[timezero; allinfo.Position(iM)];
                [token,rem] = strtok(rem,',');
                allinfo.Points(iM) = str2num(token);
                [token,rem] = strtok(rem,',');
                allinfo.Channel(iM) = str2num(token);
                allinfo.Date{iM} = '';
                
            end
        end
    else
        line=1;
    end
end

fclose(fid);    

