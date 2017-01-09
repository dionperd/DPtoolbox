%%The One Script Version 0.1
%% From Vision Analyzer Data to Coupling Data
%%Here first applied to Auditory Oddball Task

function RES=KHvaTF2CM_grid(filenum,coupletype,NTr)
%%filestr is a dummy number for the fie that has to be analysed
%%coupletype is a number:   1 for phph-coupling
%%                          2 for ph2amp-coupling
timestart=clock;

if isstr(filenum)
    filenum=str2num(filenum)
end
if isstr(coupletype)
    coupletype=str2num(coupletype)
end
if isstr(NTr)
    NTr=str2num(NTr)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%IMPORTANT GENERAL SETTINGS
%%Inputs:
if isunix
    InputPath = '/home/mpib/hosang/K_ExportOA/MlOut'; %mac
end
if ismac
    InputPath= '/Users/hosang/tardis/K_ExportOA/MlOut';
    %InputPath= '/Volumes/gridmaster2012/K_ExportOA/MlOut';
end
OutputPath=[InputPath,'/STD-VarTrial/']
if ismac || ispc
codePath = {...
    '/Users/hosang/tardis/EEG_COUPLE/io';
    %'/Volumes/InterBrain/EEGlab_VM/Denis/Software/DPtoolbox/CouplingToolbox';
    '/Users/hosang/tardis/EEG_COUPLE/CouplingToolbox';
    '/Volumes/InterBrain/EEGlab_VM/Denis/Software/DPtoolbox/DataHandling';
    '/Volumes/InterBrain/EEGlab_VM/Denis/Software/DPtoolbox/Utilities'};
end

%%the root folder name:
winRootPath = '\\MPIB10';
macRootPath = '/Volumes';

channels=[57;59;58;41;49;43;45;47;22;24;26;28;30;6;8;10;12;14;1;2;3];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%for reasons of universability:
%%%fit pathnames to OS
if ismac || ispc
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isequal(exist(OutputPath, 'dir'),7) % 7 = directory.
    mkdir(OutputPath);
end

%%so, lets change into Data Direcotory
cd(InputPath)
%cd(K_Export/MlOut)

%%Collect Data Pathes...
ii=0;
FolderDir=dir;
files=cell(2,1);
for iFile = 3:length(dir);
    if length(FolderDir(iFile).name)>=14
    if strcmpi(FolderDir(iFile).name(end-13:end),'STD-C_Coef.mat');%{'eeg','dat','bin'}
        ii=ii+1;
        files{ii}=FolderDir(iFile).name;
    end
    end
end

%%% ...and run the CM-calculation for Input-specified files
for iFile = filenum;
    disp('-----------------')
    disp(['TORTURE FILE ',files{iFile}])
    
    disp('Convert Data...')
    load(files{iFile});
    
    trials=randperm(size(Res.Coef,4),NTr);
    data=Res.Coef(:,:,:,trials);
    hdr=Res.hdr;
    
    
    %%%%Input Parameters for Coupling Toolbox
    MethodNames={'trial','trialTime','ensemble'};
    animation =[5,5,0,0];%150
    CFs=[2,4,6,8,10,12,20,30,40,50,60];%11
    PhAmpCFs=[6,48;6,60;10,60;10,50;50,10;60,10;60,6;48,6];
    PAl=[6,10];
    PAh=[48,50,60];
    PAs=[6,10,48,50,60];
    Freqs=(2:2:60);
    disp('set 21-channel system');
    if hdr.NumberOfChannels==21
        channels='all';
    else
        channels=[57;59;58;41;49;43;45;47;22;24;26;28;30;6;8;10;12;14;1;2;3];
        data=data(:,:,channels,:);
        hdr.NumberOfChannels=21;
    end
    
    %%%%%manipulate filename
    [~,fname,~]=fileparts(files{iFile});
    if fname(1:2)=='13'
        fname=['NC_',fname];
    end
    resstr=['RES-T',num2str(NTr),fname];
    resstr=strrep(strrep(strrep(strrep(resstr,'mdo',''),'1_Segm',''),'_Coef',''),'1_SegC-','');
    
    %%%set CFG for CM calculation
    cfg.fs=500;
    cfg.time=(1'/cfg.fs:1/cfg.fs:1.5)';
    %cfg.timeCalc=t;
    %cfg.Nbins = 15;
    cfg.cutTails=[0 0];
    cfg.winLen = 1/cfg.fs;
    cfg.measures={'PLV','PPC','SE','CP'};
    cfg.NstepsPwin = 1;
    cfg.stats=[];
    cfg.method=MethodNames{3};
    cfg.domain='phase';
    
    if coupletype==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% do ph2ph-CFC calculation %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for iCF1=1:length(CFs)
            for iCF2=1:1:length(CFs)
                
                cfg.fc=[CFs(iCF1),CFs(iCF2)];
                cfg.nm=[CFs(iCF1),CFs(iCF2)]/gcd(CFs(iCF1),CFs(iCF2));
                
                for iCh1=1:hdr.NumberOfChannels
                    for iCh2=1:hdr.NumberOfChannels
                        
                        %%%data input for CFC coupling
                        x1=squeeze(angle(data(:,Freqs==CFs(iCF1),iCh1,:)));
                        x2=squeeze(angle(data(:,Freqs==CFs(iCF2),iCh2,:)));
                        xIn=zeros(750,2,NTr);
                        xIn(:,1,:)=x1;
                        xIn(:,2,:)=x2;
                        
                        %%%%%%---------CALCULATE COUPLING MEASURE----%%%
                        %%%%prepare CM
                        disp(['Calculate ',strrep(num2str(cfg.fc),'  ','x'),' coupling between channel ',num2str(iCh1),' and ',num2str(iCh2)])
                        [phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds] = DPphaseCouplingPrepare(xIn,cfg);
                        %%%%calc CM
                        [PC, Dphi, cfg, statsRes] = DPphaseCoupling(phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds);
                        
                        for iMeas=1:length(cfg.measures)
                            RES.(cfg.measures{iMeas})(:,iCh1,iCh2,iCF1,iCF2)=PC{cfg.nm(1),cfg.nm(2)}.(cfg.measures{iMeas});
                            %RESULTS.(['CFC_',strrep(num2str(CFC),'  ','x')]).(cfg.measures{iMeas}){iCh1,iCh2}=PC{cfg.nm(1),cfg.nm(2)}.(cfg.measures{iMeas});
                        end
                        clear xIn;
                    end
                end
            end
        end
        save([OutputPath,resstr,'PhPh.mat'],'RES','hdr','-v7.3');
    end
    
    if coupletype==2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% do ph2amp-CFC calculation %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for iCF1=1:5
            for iCF2=1:5
                
                %%check whether numbers are multiples of each other
                if ((PAs(iCF1)*round(PAs(iCF2)/PAs(iCF1))==PAs(iCF2)) || (PAs(iCF2)*round(PAs(iCF1)/PAs(iCF2))==PAs(iCF1))) 
                    
                    cfg.fc=[PAs(min(iCF1,iCF2)),PAs(min(iCF1,iCF2))];
                    cfg.nm=[1,1];
                    
                    for iCh1=1:hdr.NumberOfChannels
                        for iCh2=1:hdr.NumberOfChannels
                            
                            %%%data input for CFC coupling
                            x1=squeeze(angle(data(:,Freqs==PAs(iCF1),iCh1,:)));
                            x2=squeeze(abs(data(:,Freqs==PAs(iCF2),iCh2,:)));
                            xIn=zeros(750,2,NTr);
                            xIn(:,1,:)=x1;
                            xIn(:,2,:)=x2;
                            
                            %%%%%%---------CALCULATE COUPLING MEASURE----%%%
                            %%%%prepare CM
                            disp(['Calculate ',strrep(num2str([PAs(iCF1),PAs(iCF2)]),'  ','x'),' ph2amp-coupling between channel ',num2str(iCh1),' and ',num2str(iCh2)])
                            [phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds] = DPphaseCouplingPrepare(xIn,cfg);
                            %%%%calc CM
                            [PC, Dphi, cfg, statsRes] = DPphaseCoupling(phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds);
                            
                            for iMeas=1:length(cfg.measures)
                                RES.(cfg.measures{iMeas})(:,iCh1,iCh2,iCF1,iCF2)=PC{cfg.nm(1),cfg.nm(2)}.(cfg.measures{iMeas});
                                %RESULTS.(['CFC_',strrep(num2str(CFC),'  ','x')]).(cfg.measures{iMeas}){iCh1,iCh2}=PC{cfg.nm(1),cfg.nm(2)}.(cfg.measures{iMeas});
                            end
                            clear xIn;
                        end
                    end
                end
            end
        end
        save([OutputPath,resstr,'PhAmp.mat'],'RES','hdr','-v7.3');
    end
            
            %save([OutputPath,resstr,'D.mat'],'data','hdr');
            disp(['Elapsed time ',num2str(etime(clock,timestart)),' s'])
            
        end