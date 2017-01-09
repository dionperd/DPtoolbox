function [h, hsub] = plotMeanSTD_2D_Net(mean,std,Measure, Colors, Groups, Conds, Ns,varargin)

Ng = max(numel(Groups),1);
Nc = max(numel(Conds),1);

%Make legend
iX=0;
if (Ng==1)
    for iC=1:Nc;
        iX=iX+1;
        xlbl{iX} = [Conds{iC}];
    end
    Ndof = Nc;
elseif (Nc==1)
    for iG=1:Ng;
        iX=iX+1;
        xlbl{iX} = [Groups{iG}];
    end
    Ndof = Ng;
else
    for iG=1:Ng;
        for iC=1:Nc;
            iX=iX+1;
            xlbl{iX} = [Groups{iG},'-',Conds{iC}];
        end
    end
    Ndof = Ng*Nc;
end

if isempty(Ns)
    titleName = [Measure,', subjects'' mean'];
    x=mean.';
else
    titleName = [Measure,', subjects'' mean / standard error'];
    %Calculate mean / standard error ratio
    for iG=1:Ng;
        for iC=1:Nc;
            x{iC,iG} = mean{iG,iC} ./ (std{iG,iC}/sqrt(Ns(iG))) ;
        end
    end
end


maxVal = -inf;
minVal = inf;
for iX=1:Ndof;
    minVal = min(minVal,min(min(x{iX})));
    maxVal = max(maxVal,max(max(x{iX})));
end

for iX=1:Ndof;
    
    if strcmpi(Measure,'ICI');
        %For ICI:
        x{iX} = DPconVec2conMatICI(x{iX},cfgMeasure);
    else
        %For PLV and CCR:
        x{iX} = DPconVec2conMat(x{iX},cfgMeasure);
    end
    
    mode = varargin{1};
    clim = varargin{2};
    if isempty(clim)
        clim = [minVal maxVal];
    end
    t = varargin{3}; %[1 2 3];
    tSel = varargin{4}; %{[1], [2], [3]};
    tlbl = varargin{5}; %{'before','during','after'};
    f = varargin{6}; %[1:10];
    fSel = varargin{7}; %{[1],[2],[3],[4],[5],[6],[7],[8],[9],[10]};
    flbl = varargin{8}; %{'2','4','6','8','10','12','14','16','18','20'};
    figTitle = varargin{9}; %{[2],'test'};
    figTitle{2} = [titleName,', ',figTitle{2}];
    
    h(iX) = figure('Name',titleName);
        
    if strcmpi(mode,'net')
        %For BrainNet connectivities:
        %BrainNetPath = '/media/dionysios/Windows/Users/dionperd/Dropbox/Dionysis/DPtoolbox/Plotting/external/BrainNetViewer/Data/'; %Denis' Linux
        %BrainNetPath = 'C:\Users\dionperd\Dropbox\Dionysis\DPtoolbox\Plotting\external\BrainNetViewer\Data'; %Denis' Windows
        surfFile = varargin{10}; %[BrainNetPafilesep,'ExampleFiles',filesep,'Custom',filesep,'EEG_21_interbrain.node'];
        optionsFile = varargin{11}; %BrainNetPath,fileseth,filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152.nv'];
        nodeFile = varargin{12}; %[BrainNetPath,p,'ExampleFiles',filesep,'Custom',filesep,'denisSymmetricInterbrainOptions.mat']; %for PLV and CCR
        %optionsFile = [BrainNetPath,filesep,'ExampleFiles',filesep,'Custom',filesep,'denisDirectedInterbrainOptions.mat'];%only for ICI
        varargin = {surfFile,nodeFile,optionsFile};
    else
        varargin = {};
    end
    [h, hsub] = plotBrainNetTimeFreq(x{iX},mode,clim,t,tSel,tlbl,f,fSel,flbl,figTitle,varargin);
    
end

        