function [EEGData, hdr, mrk, mrkPtype] = DPreadBVmatlabTRNSFRM(channels,trials,hdrFile,mrkFile)

%load workspace
T = evalin('base','whos');
for ii = 1:length(T)
   C_ =  evalin('base',[T(ii).name ';']);
   eval([T(ii).name,'=C_;']);
end
clear T C_ ii;


Nchan = size(EEGData,2); %number of channels
Ntrials = size(EEGData,3);


allChannels = 1:Nchan;
if nargin>1
    
    %Check channels
    if strcmpi(channels, 'all')
        chanindx=allChannels;
    else
        chanindx = unique(channels);
        EEGData = EEGData(:,chanindx,:);
    end
else
    chanindx=allChannels;
end
Nchan = length(chanindx);

allTrials = 1:Ntrials;
%if data is segmented into trials:
if Ntrials>1
    if nargin>2
        %Check trials
        if strcmpi(trials, 'all')
            trials=allTrials;
            
        else
            trials = unique(trials);
            EEGData = EEGData(:,:,trials);
        end
    else
        trials=allTrials;
    end
end
Ntrials = length(Ntrials);



if nargin<3
    hdrFile = '';
elseif ~exist(hdrFile,'file')
    hdrFile='';
end

if isempty(hdrFile) %if we read hdr from MATLAB workspace
    
    
    hdr.DataFile = [FileName,'.eeg'];
    hdr.MarkerFile = [FileName,'vmrk'];
    hdr.DataFormat = 'BINARY';
    hdr.DataOrientation = 'MULTIPLEXED';
    hdr.DataType = 'TIMEDOMAIN';
    hdr.NumberOfChannels = Nchan;
    hdr.DataPoints = Ntrials*length(EEGTime);
    hdr.SamplingInterval = Properties.SamplingInterval;
    hdr.Fs = 10^6/hdr.SamplingInterval; %sampling frequency
    hdr.time = EEGTime;
    
    if isfield(Properties, 'SegmentationType')
        if Properties.SegmentationType==1
            hdr.SegmentationType = 'MARKERBASED';
        end
    end
    if isfield(Properties, 'SegmentDataPoints')
        hdr.SegmentDataPoints = Properties.SegmentSize;
    end
    
    hdr.BinaryFormat = 'IEEE_FLOAT_32';
    
    iC=0;
    for ii = chanindx;
        iC=iC+1;
        hdr.label{iC} = Properties.Channels(ii).Name;
        hdr.reference{iC} = Properties.Channels(ii).RefName;
        hdr.units{iC} ='?V';
        hdr.coords.radius(iC) = Properties.Channels(ii).CoordsRadius;
        hdr.coords.theta(iC) = Properties.Channels(ii).CoordsTheta;
        hdr.coords.phi(iC) = Properties.Channels(ii).CoordsPhi;
    end
    
else %if we read hdr from a .vhdr file
    
    hdr=DPreadBVvhdr(hdrFile);
    hdr.time=EEGTime;
    
    %Reduce hdr if a specific selection of channels and/or trials has been
    %made:
    if ~isequal(chanindx,allChannels)
        hdr.NumberOfChannels = size(data,2);
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
        hdr.nTrials = size(data,3);
        hdr.DataPoints = size(data,1)*hdr.nTrials ;
        hdr.nSamples = hdr.DataPoints;
    end
end


%Read markers
if nargin<4
    mrkFile = '';
elseif ~exist(mrkFile,'file')
    mrkFile='';
end

if isempty(mrkFile) %if we read mrk from MATLAB workspace
    
    MrkTypes = {'New Segment','Time 0','Bad Interval','Stimulus','Response','Marker','Comment'};
    Nmrks=length(Markers); %number of markers
    
    for iM=1:Nmrks;
        mrk.type{iM} = Markers(iM).Type;
        if any(strcmpi(mrk.type{iM},MrkTypes))
            mrk.label{iM} = Markers(iM).Description;
        else
            mrk.label{iM} = '';
        end
        mrk.Position(iM) = Markers(iM).Position+1;
        mrk.Points(iM)=Markers(iM).Points;
        mrk.Channel(iM)=Markers(iM).Channel;
        if (mrk.Channel(iM)==-1)
            mrk.Channel(iM)=0;
        end
        if strcmpi(mrk.type{iM},'New Segment')
            mrk.Date{iM} = [Markers(iM).Date(7:10),Markers(iM).Date([4,5]),Markers(iM).Date([1,2]),Markers(iM).Date([12,13]),Markers(iM).Date([15,16]),Markers(iM).Date([18,19]),'000000'];
        else
            mrk.Date{iM} ='';
        end
        
    end
    
else %if we read mrk from a .vmrk file
    
     mrk = DPreadBVvmrk(mrkFile);

end


%Reduce mrk if a specific selection of trials has been
%made:
if ~isequal(trials,allTrials)
    
    Nmrks = length(mrk.type); %number of markers in the old structure
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
        
        %check if this trial is included in the required set
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
            mrkPtype.MrkTypSTR{iT}.Date{ind} = mrk.Date{iM};
        end
        
    end
    
end


