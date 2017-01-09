function DPrunPLS_PSD_JLDmanyAnal

Analyses = {'FCD'};%,'JLD','COS','JLDn'};...;%{'FCD','FCDlogit','JLD','JLDlog',...
   %'FCDsurr','FCDlogitSurr','JLDsurr','JLDlogSurr'};

Nanal = length(Analyses);

MainPath = 'E:\home\dionysios\Documents\PostDocMRS\Results\EEG_PSD_JLD\Results_EEG_PSDfinal_FCD';
%MainPath = '/Users/dionysos/Dropbox/Dionysis/PostDocMRS/Results/EEG_PSD_JLD/';
%MainPath ='C:\Users\CoordAgeEEG\Dropbox\Results\EEG_PSD_JLD\PSDnorm';

for iN = 1:Nanal;
    
%     thisPath = [MainPath,filesep,Analyses{iN}];
    thisPath = MainPath;
    thisPLSPath = [thisPath,filesep,'StatsLifespanTask'];
    
    if ~exist(thisPLSPath,'dir')
        mkdir(thisPLSPath)
    end
    
    DPmakePLSmatPSD_JLDfun(thisPath,thisPLSPath);
    
   %DPmakePLSmatPSD_JLDfun_aging(thisPath,thisPLSPath);

    DPrunPLS_PSD_JLDnewFun(thisPLSPath);
 
    %movefile([thisPLSPath,filesep,'PLSres.fig'],[thisPLSPath,filesep,'PLSres1.fig'])
    %DPplotPLS_PSD_JLDtogetherFun(thisPLSPath);
    DPplotPLS_PSD_JLDfun(thisPLSPath);
    
    %DPplotPLS_PSD_JLDfun_aging(thisPLSPath);

end