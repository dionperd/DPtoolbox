function DPrunPLS_PSDnorm_manyAnal

%Analyses = {'FCD','COS','JLD','JLDn'};...;%{'FCD','FCDlogit','JLD','JLDlog',...
   %'FCDsurr','FCDlogitSurr','JLDsurr','JLDlogSurr'};

%Nanal = length(Analyses);

MainPath1 ='E:\home\dionysios\Documents\PostDocMRS\Results\EEG_PSD_JLD\PSnorm';
%MainPath = 'C:\Users\dionperd\Dropbox\Dionysis\PostDocMRS\Results\EEG_PSD_JLD';
%MainPath = '/Users/dionysos/Dropbox/Dionysis/PostDocMRS/Results/EEG_PSD_JLD';
MainPath2 ='E:\home\dionysios\Documents\PostDocMRS\Results\EEG_PSD_JLD\PSnorm';

%for iN = 1:Nanal;
    
    %thisPath = [MainPath,filesep,Analyses{iN}];
    thisPLSPath = [MainPath2,filesep,'StatsAging'];
    
    if ~exist(thisPLSPath,'dir')
        mkdir(thisPLSPath)
    end
    
%      DPmakePLSmatPSDnormfun(MainPath1,thisPLSPath);
DPmakePLSmatPSDnormfun_aging(MainPath1,thisPLSPath);
% %     
      DPrunPLS_PSDnormFun(thisPLSPath);
    
    %DPplotPLS_PSD_JLDtogetherFun(thisPLSPath);
   % DPplotPLS_PSDnormfun(thisPLSPath);
%end