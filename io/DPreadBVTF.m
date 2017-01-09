function [data hdr mrk stimulusMRK, responseMRK, segmentMRK, timezeroMRK] = DPreadBVTF(filename,channels,trials)


[pathstr, name, ~] = fileparts(filename); 
%older version
% %get the position of the last '.' in the filename
% pointInd = strfind(filename,'.');
% pointInd=pointInd(end);
% %this is the filename without the extension
% filenameNoExt = filename(1:pointInd);

hdrFilename = fullfile(pathstr,[name,'.vhdr']);
if exist(hdrFilename,'file')
    
    %Read header
    hdr = DPreadBVvhdrTF(hdrFilename);
    allChannels = 1:hdr.NumberOfChannels;
    
    if nargin>1
        
        %Check channels
        if strcmpi(channels, 'all')
            chanindx=allChannels;
        else
            chanindx = unique(channels);
        end
    else
        chanindx=allChannels;
    end
    
    
    allTrials = 1:hdr.nTrials;
    %if data is segmented into trials:
    if hdr.nTrials>1
        
        if nargin>2
            %Check trials
            if strcmpi(trials, 'all')
                trials=allTrials;
                
            else
                trials = unique(trials);
            end
        else
            trials=allTrials;
        end
    else
        trials=1;
    end
    
    %Read data
    data = DPreadBVeegTF(filename, hdr,chanindx,trials);
    
else
    error('there is no header file in the specified path')
end

%Reduce hdr if a specific selection of channels and/or trials has been
%made:
if ~isequal(chanindx,allChannels)
    hdr.NumberOfChannels = size(data,3);
    hdr.label = hdr.label(chanindx);
    hdr.reference = hdr.reference(chanindx);
    if isfield(hdr,'units')
        hdr.units = hdr.units(chanindx);
    end
    if isfield(hdr,'resolution')
        hdr.resolution = hdr.resolution(chanindx);
    end
    hdr.coords.radius = hdr.coords.radius(chanindx);
    hdr.coords.theta = hdr.coords.theta(chanindx);
    hdr.coords.phi = hdr.coords.phi(chanindx);
end
if ~isequal(trials,allTrials)
    hdr.nTrials = size(data,4);
    hdr.DataPoints = size(data,1)*hdr.nTrials ;
    hdr.nSamples = hdr.DataPoints;  
end


%Read markers
mrkFilename = fullfile(pathstr,[name,'.vmrk']);
if exist(mrkFilename,'file')
    
    mrk = DPreadBVvmrk(mrkFilename);
    Nmrks = length(mrk.type); %number of markers in the old structure

    %Reduce mrk if a specific selection of trials has been
    %made:
    if ~isequal(trials,allTrials)
        
        iT=0; %meter of New Segment markers
        %Find the indices of all New Segment, i.e. trials', markers
        for iM=1:Nmrks;
            if strcmpi(mrk.type{iM},'New Segment')
                iT=iT+1;
                newSegmINDs(iT) =iM;
            end
        end
        
        Nsegms = length(newSegmINDs); %number of new segments/trials
        oldmrk = mrk;%store the old mrk
        mrk=struct();%initialize the new one
        NmrksNew = 0;%the initial number of markers in the new structure
        
        %Select only the markers of the included trials
        for iT=1:Nsegms-1;
            
            %check of this trial is included in the required set
            if ismember(iT, trials)
                
                %if it is not the last one
                if iT<Nsegms
                    %choose all markers until the next new segment
                    thisMrks = newSegmINDs(iT+1) - newSegmINDs(iT)-1;
                else%if it is the last one
                    %choose all markers until the last one
                     thisMrks = Nmrks - newSegmINDs(iT);
                end
                
                %add to the new mrk structure
                mrk.type(NmrksNew+1:NmrksNew+thisMrks+1) = oldmrk.type(newSegmINDs(iT):newSegmINDs(iT)+thisMrks);
                mrk.label(NmrksNew+1:NmrksNew+thisMrks+1) = oldmrk.label(newSegmINDs(iT):newSegmINDs(iT)+thisMrks);
                mrk.Position(NmrksNew+1:NmrksNew+thisMrks+1) = oldmrk.Position(newSegmINDs(iT):newSegmINDs(iT)+thisMrks);
                mrk.Points(NmrksNew+1:NmrksNew+thisMrks+1) = oldmrk.Points(newSegmINDs(iT):newSegmINDs(iT)+thisMrks);
                mrk.Channel(NmrksNew+1:NmrksNew+thisMrks+1) = oldmrk.Channel(newSegmINDs(iT):newSegmINDs(iT)+thisMrks);
                mrk.Date(NmrksNew+1:NmrksNew+thisMrks+1) = oldmrk.Date(newSegmINDs(iT):newSegmINDs(iT)+thisMrks);
                NmrksNew = length(mrk.type);
            end
            
        end
        Nmrks=NmrksNew;
 
    end
    
else
    mrk = [];
end


if (nargout>3)
    
    MrkTypes = {'New Segment','Time 0','Bad Interval','Stimulus','Response','Marker','Comment'};
    MrkTypSTR = {'NewSegment','Time0','BadInterval','Stimulus','Response','Marker','Comment'};
    Nt=length(MrkTypes);
    MrkTypIND = MrkTypSTR; %initialization
    
    for iT = 1:Nt;
        MrkTypIND{iT} = ['i',MrkTypSTR{iT}];
        mrkPtype.(MrkTypIND{iT}) = 0;
        mrkPtype.(MrkTypSTR{iT}) = [];
    end
    
    if (~isempty(mrk))
        
        %reconstruct according to the new mrk
        for iM=1:Nmrks;
            
            for iT = 1:Nt;
                if strcmpi(mrk.type{iM},MrkTypes{iT})
                    mrkPtype.(MrkTypIND{iT})=mrkPtype.(MrkTypIND{iT})+1;
                    ind=mrkPtype.(MrkTypIND{iT});
                    break;
                end
            end
            mrkPtype.MrkTypSTR{iT}.label{ind} = mrk.label{iM};
            mrkPtype.MrkTypSTR{iT}.Position(ind) = (mrk.Position(iM)-1)/hdr.Fs*1000;
            mrkPtype.MrkTypSTR{iT}.Points(ind) = mrk.Points(iM);
            mrkPtype.MrkTypSTR{iT}.Channel(ind) = mrk.Channel(iM);
            %mrkPtype.MrkTypSTR{iT}.Date{ind} = mrk.Date{iM}; %Only New Segment have Date field...
        end
        
    end
    
end


