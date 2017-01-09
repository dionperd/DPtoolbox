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
%     if isfield(hdr,'resolution')
%         hdr.resolution = hdr.resolution(chanindx);
%     end
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
    
    %Finally reconstruct stimulusMRK, responseMRK, segmentMRK & timezeroMRK
    %initialize
    iSeg = 0;
    iStim = 0;
    iResp =0;
    iTim0 =0;
    segmentMRK = [];
    stimulusMRK={};
    responseMRK={};
    timezeroMRK=[];
    %reconstruct according to the new mrk
    for iM=1:Nmrks;
        if strcmpi(mrk.type{iM},'New Segment')
            iSeg=iSeg+1;
            segmentMRK(iSeg)=(mrk.Position(iM)-1)/hdr.Fs*1000;
        elseif strcmpi(mrk.type{iM},'Stimulus')
            iStim=iStim+1;
            stimulusMRK{iStim,1} = str2num(mrk.label{iM}(2:end));
            stimulusMRK{iStim,2} = (mrk.Position(iM)-1)/hdr.Fs*1000;
        elseif strcmpi(mrk.type{iM},'Response')
            iResp=iResp+1;
            responseMRK{iResp,1} = str2num(mrk.label{iM}(2:end));
            responseMRK{iResp,2} = (mrk.Position(iM)-1)/hdr.Fs*1000;
        elseif  strcmpi(mrk.type{iM},'Time 0')
            iTim0=iTim0+1;
            timezeroMRK(iTim0)=(mrk.Position(iM)-1)/hdr.Fs*1000;
        end
    end
    

else
    
    segmentMRK = [];
    stimulusMRK={};
    responseMRK={};
    timezeroMRK=[];
end



