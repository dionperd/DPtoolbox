function [ output_args ] = KHcoupleEEG( input_args )
%%Take Time-Frequency Data from Vision Analyzer and calculate intra-subject CFC


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%IMPORTANT GENERAL SETTINGS
%%Inputs:
InputPath = '/Volumes/LIP-Git2/Denis/Karl/K_Export/MlOut'; %mac
OutputPath=[InputPath,'\MlOut']
codePath = {...
            '/Volumes/InterBrain/EEGlab_VM/Denis/Software/DPtoolbox/io';
            '/Volumes/InterBrain/EEGlab_VM/Denis/Software/DPtoolbox/CouplingToolbox'};

%%the root folder name:
winRootPath = '\\MPIB10';
macRootPath = '/Volumes';

channels=[57;59;58;41;49;43;45;47;22;24;26;28;30;6;8;10;12;14;1;2;3];
trials='all';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%for reasons of universability:
%%%fit pathnames to OS
for ii=1:length(codePath);
    if ispc
        codePath{ii} = strrep(codePath{ii}, macRootPath, winRootPath);
        InputPath = strrep(InputPath, macRootPath, winRootPath);
        InputPath = strrep(InputPath, '/', '\');
        codePath{ii} = strrep(codePath{ii}, '/', '\');
        OutputPath = strrep(OutputPath, macRootPath, winRootPath);
        OutputPath = strrep(OutputPath, '/', '\');
    end
    if ismac
        codePath{ii} = strrep(codePath{ii}, winRootPath,macRootPath);
        codePath{ii} = strrep(codePath{ii}, '\', '/');
        OutputPath = strrep(OutputPath, winRootPath, macRootPath);
        OutputPath = strrep(OutputPath, '\', '/');
        InputPath = strrep(InputPath, winRootPath, macRootPath);
        InputPath = strrep(InputPath, '\', '/');
    end
    if isunix
        codePath{ii} = strrep(codePath{ii}, '\', '/');
    end
    addpath(codePath{ii})
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Input Parameters for Coupling Toolbox
MethodNames={'trial','trialTime','ensemble'};
animation =[5,5,0,0];%150
CFCs=[2,10;6,60;6,48;10,50;10,60;10,40];


%%%set CFG for CM calculation
        cfg.fs=500;
        cfg.time=(0:1/cfg.fs:1.5);
        %cfg.timeCalc=t;
        %cfg.Nbins = 15;
        cfg.cutTails=[0 0];
        %cfg.winlen = 1/fs;
        cfg.measures={'PLV','IC','PPC','PLI','SE','MI','CP'};
        cfg.NstepsPwin = 1;
        cfg.stats=[];
        cfg.method=MethodNames{1};
        cfg.nm=[1 5];

%%so, lets change into Data Direcotory
cd(InputPath)
        
%%Collect Data Pathes
ii=0;
FolderDir=dir;
files=cell(2,1);
for iFile = 3:length(dir);
    if any(strcmpi(FolderDir(iFile).name(end-2:end),{'eeg','dat','bin'}));%{'eeg','dat','bin'}
        ii=ii+1;
        files{ii}=FolderDir(iFile).name;
        disp('-----------------')
        disp(['TORTURE FILE ',FolderDir(iFile).name])
%         KHVAcouplefun(cfg,channels,trials,FolderDir(iFile).name,winRootPath,macRootPath);
    end
end


function KHVAcouplefun(procfun,cfg,channels,trials,filename,winRootPath,macRootPath)
%%% Read Vision Analyzer Data
[data hdr mrk stimulusMRK, responseMRK, segmentMRK, timezeroMRK] = KHreadBVTF(filename,channels,trials)

%%%data input for coupling toolbox
xIn=[];



%%%%%%---------CALCULATE COUPLING MEASURE----%%%
%%%%prepare CM
[phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds] = DPphaseCouplingPrepare(xIn,cfg);
%%%%calc CM
[PC, Dphi, cfg, statsRes] = DPphaseCoupling(phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds);
