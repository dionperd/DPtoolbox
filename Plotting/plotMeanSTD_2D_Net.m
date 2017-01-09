function [h, hsub] = plotMeanSTD_2D_Net(mean,std,Measure, Colors, Groups, Conds, Ns,h,varargin)

Ng = max(numel(Groups),1);
Nc = max(numel(Conds),1);

cfgMeasure = varargin{end-1};
type = varargin{end};

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



for iX=1:Ndof;
    
    mode = varargin{1};
    clim = varargin{2};
    t = varargin{3}; %[1 2 3];
    tSel = varargin{4}; %{[1], [2], [3]};
    tlbl = varargin{5}; %{'before','during','after'};
    f = varargin{6}; %[1:10];
    fSel = varargin{7}; %{[1],[2],[3],[4],[5],[6],[7],[8],[9],[10]};
    flbl = varargin{8}; %{'2','4','6','8','10','12','14','16','18','20'};
    figTitle = varargin{9}; %{[2],'test'};
    %figTitle{2} = [titleName,', ',figTitle{2}];
    
    if strcmpi(Measure,'ICI');
        %For ICI:
        if strcmpi(type,'within')
            if isempty(clim)
                [x{iX}, clim] = DPconVec2conMatICIwithin(x{iX},cfgMeasure);
            else
                x{iX} = DPconVec2conMatICIwithin(x{iX},cfgMeasure);
            end
        elseif strcmpi(type,'between')
            if isempty(clim)
                [x{iX}, clim] = DPconVec2conMatICIbet(x{iX},cfgMeasure);
            else
                x{iX} = DPconVec2conMatICIbet(x{iX},cfgMeasure);
            end
            else
            if isempty(clim)
                [x{iX}, clim] = DPconVec2conMatICI(x{iX},cfgMeasure);
            else
                x{iX} = DPconVec2conMatICI(x{iX},cfgMeasure);
            end
        end
    else
        %For PLV and CCR:
        if strcmpi(type,'within')
            if isempty(clim)
                [x{iX}, clim] = DPconVec2conMatWithin(x{iX},cfgMeasure);
            else
                x{iX} = DPconVec2conMatWithin(x{iX},cfgMeasure);
            end
        elseif strcmpi(type,'between')
            if isempty(clim)
                [x{iX}, clim] = DPconVec2conMatBet(x{iX},cfgMeasure);
            else
                x{iX} = DPconVec2conMatBet(x{iX},cfgMeasure);
            end
            else
            if isempty(clim)
                [x{iX}, clim] = DPconVec2conMat(x{iX},cfgMeasure);
            else
                x{iX} = DPconVec2conMat(x{iX},cfgMeasure);
            end
        end
    end
    
    
    
    if strcmpi(mode,'net')
        
        figTitle{2} = [titleName,', ',figTitle{2},' ', xlbl{iX}];
        
        %For BrainNet connectivities:
        %BrainNetPath = '/media/dionysios/Windows/Users/dionperd/Dropbox/Dionysis/DPtoolbox/Plotting/external/BrainNetViewer/Data/'; %Denis' Linux
        %BrainNetPath = 'C:\Users\dionperd\Dropbox\Dionysis\DPtoolbox\Plotting\external\BrainNetViewer\Data'; %Denis' Windows
        surfFile = varargin{10}; %[BrainNetPafilesep,'ExampleFiles',filesep,'Custom',filesep,'EEG_21_interbrain.node'];
        optionsFile = varargin{11}; %BrainNetPath,fileseth,filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152.nv'];
        nodeFile = varargin{12}; %[BrainNetPath,p,'ExampleFiles',filesep,'Custom',filesep,'denisSymmetricInterbrainOptions.mat']; %for PLV and CCR
        %optionsFile = [BrainNetPath,filesep,'ExampleFiles',filesep,'Custom',filesep,'denisDirectedInterbrainOptions.mat'];%only for ICI
        [h(iX), hsub] = plotBrainNetTimeFreq_Net(0,x{iX},type,mode,clim,t,tSel,tlbl,f,fSel,flbl,figTitle,surfFile,nodeFile,optionsFile);
        set(h(iX),'Name',titleName);
    else
        figTitle{2} = [titleName,', ',figTitle{2},' ', xlbl{iX}];
        
        if strcmpi(Measure,'ICI')
            [h(iX), hsub] = plotBrainNetTimeFreq_MapICI(0,x{iX},type,clim,t,tSel,tlbl,f,fSel,flbl,figTitle);
        else
            [h(iX), hsub] = plotBrainNetTimeFreq_Map(0,x{iX},type,clim,t,tSel,tlbl,f,fSel,flbl,figTitle);
        end
        set(h(iX),'Name',titleName);
    end
    
end

        