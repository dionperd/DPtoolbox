function plotBrainLVs_2D_Net(LV,BS,Measure,cfgMeasure,h,opts)

if isempty(BS)
    x = LV;
elseif strcmpi(opts.method,'BootsRatios')
    x = BS;
else
    x = LV;
end
    
mode = varargin{1};
clim = varargin{2};
t = varargin{3}; %[1 2 3]; 
tSel = varargin{4}; %{[1], [2], [3]};
tlbl = varargin{5}; %{'before','during','after'};
f = varargin{6}; %[1:10];
fSel = varargin{7}; %{[1],[2],[3],[4],[5],[6],[7],[8],[9],[10]};
flbl = varargin{8}; %{'2','4','6','8','10','12','14','16','18','20'};
figTitle = varargin{9}; %{[2],'test'};
type = varargin{end};

if strcmpi(Measure,'ICI');
    %For ICI:
    if strcmpi(type,'within')
        if isempty(clim)
            [x, clim] = DPconVec2conMatICIwithin(x,cfgMeasure);
        else
            x = DPconVec2conMatICIwithin(x,cfgMeasure);
        end
    elseif strcmpi(type,'between')
        if isempty(clim)
            [x, clim] = DPconVec2conMatICIbet(x,cfgMeasure);
        else
            x = DPconVec2conMatICIbet(x,cfgMeasure);
        end
    else
        if isempty(clim)
            [x, clim] = DPconVec2conMatICI(x,cfgMeasure);
        else
            x = DPconVec2conMatICI(x,cfgMeasure);
        end
    end
else
    %For PLV and CCR:
    if strcmpi(type,'within')
        if isempty(clim)
            [x, clim] = DPconVec2conMatWithin(x,cfgMeasure);
        else
            x = DPconVec2conMatWithin(x,cfgMeasure);
        end
    elseif strcmpi(type,'between')
        if isempty(clim)
            [x, clim] = DPconVec2conMatBet(x,cfgMeasure);
        else
            x = DPconVec2conMatBet(x,cfgMeasure);
        end
    else
        if isempty(clim)
            [x, clim] = DPconVec2conMat(x,cfgMeasure);
        else
            x = DPconVec2conMat(x,cfgMeasure);
        end
    end
end



if strcmpi(mode,'net')
    %For BrainNet connectivities:
    %BrainNetPath = '/media/dionysios/Windows/Users/dionperd/Dropbox/Dionysis/DPtoolbox/Plotting/external/BrainNetViewer/Data/'; %Denis' Linux
    %BrainNetPath = 'C:\Users\dionperd\Dropbox\Dionysis\DPtoolbox\Plotting\external\BrainNetViewer\Data'; %Denis' Windows
    surfFile = varargin{10}; %[BrainNetPafilesep,'ExampleFiles',filesep,'Custom',filesep,'EEG_21_interbrain.node'];
    optionsFile = varargin{11}; %BrainNetPath,fileseth,filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152.nv'];
    nodeFile = varargin{12}; %[BrainNetPath,p,'ExampleFiles',filesep,'Custom',filesep,'denisSymmetricInterbrainOptions.mat']; %for PLV and CCR
    %optionsFile = [BrainNetPath,filesep,'ExampleFiles',filesep,'Custom',filesep,'denisDirectedInterbrainOptions.mat'];%only for ICI
    [h, hsub] = plotBrainNetTimeFreq_Net(h,x,type,mode,clim,t,tSel,tlbl,f,fSel,flbl,surfFile,nodeFile,optionsFile);
else
    if strcmpi(Measure,'ICI');
        [h, hsub] = plotBrainNetTimeFreq_MapICI(h,x,type,clim,t,tSel,tlbl,f,fSel,flbl,figTitle);
    else
        [h, hsub] = plotBrainNetTimeFreq_Map(h,x,type,clim,t,tSel,tlbl,f,fSel,flbl,figTitle);
    end
end

