function DPmakePLSmatPSDnormfun_aging(RESfolder,PLSfolder)
%function [PLSres, PLScfg] = DPrunPLStask( PLSfolder, RESfolder,  RESstruct, Conds, Groups, Measures, MeasOutsStruct, Ntr, Ns, Subjs)
%DPrunPLS This function calculates PLS analysis for several measures on a
%potentially multigroup and multisubject dataset

%   Outputs:
% -PLSres: structure with PLS results for each measure
% -PLScfg: structure with PLS configuration (PLSmat included)  for each measure

%% %---------------------A. Stack PLS matrix----------------------------------
%   Inputs:
% RESfolder='./Results';
% %-RESfolder: full path to the folder of results .mat files
% PLSmatFolder = 'C:\Users\dionperd\Dropbox\Dionysis\PostDocMRS\Results\EEG_PSD_JLD\PLSallParamsTogether';
PLSmatFolder = RESfolder;
StatsFolder = PLSfolder;%[PLSmatFolder,filesep,'Stats'];
if ~exist(StatsFolder,'dir')
    mkdir(StatsFolder)
end
PLSmatFolder = StatsFolder;
%-PLSmatFolder: full path to the folder where PLS matrix will be saved
RESstruct = 'PSDstats';
%-RESstruct: name of the result structure (string)
Conds={'RestEC','Odd_nC','Odd_C'};%'RestEC','Odd_nC','Odd_C'
%-Conds: cellstring of the names of conditions
Groups={'mdo13','mdo14'};% 'mdo11','mdo12',
%-Groups: cellstring of the names of Groups 
%(a string common to all files, e.g., the extension, if there is only 1 group)
Ns= [31 28]; %   % RestEO 24 28 
%Ns= [30 26]; %  23 28 % The rest conditions without outliers
%-Ns: number of subjects of each group (vector of positive integers)
%load('C:\Users\CoordAgeEEG\Dropbox\Results\Subjects.mat')
Subjs = {};%
%-Subj: cell of cellstrings of the names of subjects for each group in case we want to use a
%       specific subset of them, otherwise Subjs={}
Measures = {'norm_m','norm_s','norm_k'};%,'norm_m','norm_s','norm_ma','norm_sa' ,'ma','sa','ka',
%-Measures: a cellstring of the names of measures to be calculated
%                   all 4 JLD params   
%                  scales>1.004sec 10 freqs 4 params  
MeasOutsStruct =  {[1 10],[1 10],[1 10]};%{[36 10 4]};[171 10],[171 10],[171 10],[171 10]

%-MeasOutsStruct: a cell of the same size as Measures with vectors
%                 corresponding to the size of the specific measure's 
%                 result PER TRIAL or for the mean data,
%                 but the trial index has to be the LAST one, if any
Ntr=1;
%-Ntr: number of trials to be considered for the within subject average (integer),
%      if Ntr = 0, an average over all trials is taken
%      if Ntr = 1, only the first trial will be taken, which is equivalent
%                  to taking mean data already calculated
%
%---------------------For Plotting-----------------------------------------
%Inputs for plotting:
PlotGroups={'YA','OA'};%'YC','OC',
PlotConds={'REC','UOT','AOT'}; %'RestEO',
%Colors = {rgb('DodgerBlue'),rgb('ForestGreen'),rgb('Red'), rgb('SteelBlue'),rgb('Teal'),rgb('DarkRed')};
Colors = {...
          %rgb('Magenta'),rgb('DarkOrchid'),rgb('Indigo'),...
          %rgb('Crimson'),rgb('IndianRed'),rgb('Maroon')...
          rgb('Lime'),rgb('ForestGreen'),rgb('DarkGreen'),...
          rgb('Cyan'),rgb('SteelBlue'),rgb('DarkBlue')...
          };
%AllMeasures = {'\mu','\sigma','\alpha','\beta'};
%varargin={}; %any other inputs you would like to give to the plotting function
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds, Ns)plotMeanSTD_2D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
%plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns)plotMeanSTD_1D(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,varargin);
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds)plotMeanSTD_0D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
%method = 'MeanStdErr';
%%plotMeasure = @(M,STD,Ns,Ng,Nc,Groups,Conds,method,measure,cfg)plotUnivComplmapMeanSTD(M,STD,Ns,Ng,Nc,PlotGroups,PlotConds,method,measure,cfg);
plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns,cfgMeasure)plotMeanSTD_2D_JLDparams(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,cfgMeasure);
%plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns,cfgMeasure)plotMeanSTD_2D_JLDparamsAlltogether(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,cfgMeasure);
%plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns,cfgMeasure)plotMeanSTD_1D_JLDparamsPerFreq(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,cfgMeasure);
% % -plotMeasure: a handle to a function to plot each measure
%---------------------For Plotting-----------------------------------------
%----------------------------End of Inputs---------------------------------




Nc = max(numel(Conds),1);
Ng = max(numel(Groups),1);
Nm =  numel(Measures);


Nr = Ns*Nc; %number of rows of PLSmat
%Initialize PLS matrices for each (sub)measure
for iM = 1:Nm;%For each measure...
    %Initialize
    PLSmat.(Measures{iM}) = cell(Ng,1);
    for iG=1:Ng;  %...for each group...
        PLSmat.(Measures{iM}){iG} = nan(Nr(iG),prod(MeasOutsStruct{iM}));
    end
    PLSmatStats.(Measures{iM}).mean = cell(Ng,Nc);
    PLSmatStats.(Measures{iM}).std = cell(Ng,Nc);
end

cd(RESfolder)

D = dir;
Nf=numel(D);
%The function to choose among the files...
if (isempty(Subjs))
    %...without subjects selection...
    ChooseFileCommand = '( ~isempty(strfind(D(iF).name,Groups{iG})) && ~isempty(strfind(D(iF).name,Conds{iC})) )';
else
    %...and with subjects selection...
    ChooseFileCommand = '( ~isempty(strfind(D(iF).name,Groups{iG})) && ~isempty(strfind(D(iF).name,Conds{iC}))&& ~isempty(strfind(D(iF).name,Subjs{iG}{iS})) )';
end

%A set of functions to get the data for each row of PLS mat:
if (Ntr==1) %this is the case where we already have 1 value for all trials (either mean or ensemble)
    for iM = 1:Nm;%...for each measure...
        %...get into temp the already averaged or ensemble value
        %for this subject...
        tempCommand{iM} = [RESstruct,'.',Measures{iM}];
    end
    
elseif (Ntr==0)%this is the case where we need to calculate the mean of ALL trials
    for iM = 1:Nm;%...for each measure...
        Mdim(iM) = length(MeasOutsStruct{iM});
        %...get into temp the average across all trials
        %for this subject...
        tempCommand{iM} = ['mean(',RESstruct,'.',Measures{iM},', Mdim(iM)+1)'];
    end
    
else
    for iM = 1:Nm;%...for each measure...
        Mdim(iM) = length(MeasOutsStruct{iM});
        %...construct the string that will index the
        %subset of selected trials...
        str{iM}='(';
        for iD = 1:Mdim(iM);
            str{iM} = [str{iM},':,'];
        end
        str=[str{iM},'[1:Ntr])'];
        %...get into temp the average across the first Ntr trials
        %for this subject...
        tempCommand{iM} = ['mean(',RESstruct,'.',Measures{iM},str{iM},', Mdim(iM)+1)'];
    end
end %if Ntr


tic
%Stack PLS matrix for each measure
disp('Stacking PLS matrix for each measure...')
for iG = 1:Ng; %For each group...
    disp('Group=')
    disp(Groups{iG})
    for iC = 1:Nc; %...for each condition...
        disp('Condition=')
        disp(Conds{iC})
        
        for iS=1:Ns(iG); %...for each subject...
            disp('Subject ')
            disp(iS)
            for iF=1:Nf; %...for each file...
                
                try
                %...if this file corresponds to this group, condition (and subject)...
                if eval(ChooseFileCommand)

                    %...load the file...
                    load(D(iF).name)
                    
                    disp('PLSmat row index: ')
                    RowIND = (iC-1)*Ns(iG)+iS;
                    disp(RowIND)
                    
                    PSDstats.norm_m = DPcheckMeasureForPLSmat(PSDstats.norm_m,-1,1,D(iF).name,Groups{iG},Conds{iC},'norm_m',iG,iC,iS,iF);
                    
                    PSDstats.norm_s = DPcheckMeasureForPLSmat(PSDstats.norm_s,0,1,D(iF).name,Groups{iG},Conds{iC},'norm_s',iG,iC,iS,iF);
                    
                    PSDstats.norm_k = DPcheckMeasureForPLSmat(PSDstats.norm_k,0,10,D(iF).name,Groups{iG},Conds{iC},'norm_k',iG,iC,iS,iF);
% 
%                     PSDstats.norm_ma = DPcheckMeasureForPLSmat(PSDstats.norm_ma,-1,1,D(iF).name,Groups{iG},Conds{iC},'norm_ma',iG,iC,iS,iF);
%                     
%                     PSDstats.norm_sa = DPcheckMeasureForPLSmat(PSDstats.norm_sa,0,1,D(iF).name,Groups{iG},Conds{iC},'norm_sa',iG,iC,iS,iF);
%                     
%                     PSDstats.norm_ka = DPcheckMeasureForPLSmat(PSDstats.norm_ka,0,10,D(iF).name,Groups{iG},Conds{iC},'norm_ka',iG,iC,iS,iF);

%                     PSDstats.m = DPcheckMeasureForPLSmat(PSDstats.m,-1,1,D(iF).name,Groups{iG},Conds{iC},'m',iG,iC,iS,iF);
%                     
%                     PSDstats.s = DPcheckMeasureForPLSmat(PSDstats.s,0,1,D(iF).name,Groups{iG},Conds{iC},'s',iG,iC,iS,iF);
% 
%                     PSDstats.k = DPcheckMeasureForPLSmat(PSDstats.k,0,10,D(iF).name,Groups{iG},Conds{iC},'k',iG,iC,iS,iF);

%                     PSDstats.ma = DPcheckMeasureForPLSmat(PSDstats.ma,-10,10,D(iF).name,Groups{iG},Conds{iC},'ma',iG,iC,iS,iF);
%                     
%                     PSDstats.sa = DPcheckMeasureForPLSmat(PSDstats.sa,0,10,D(iF).name,Groups{iG},Conds{iC},'sa',iG,iC,iS,iF);
% 
%                     PSDstats.ka = DPcheckMeasureForPLSmat(PSDstats.ka,0,10,D(iF).name,Groups{iG},Conds{iC},'ka',iG,iC,iS,iF);


                    
                    for iM = 1:Nm;%...for each measure...
                        
                        %...get this subject's data into temp...
                        temp=eval(tempCommand{iM});
                        
                        %                         disp(Measures{iM})
                        %                         disp(size(temp))
                        
                        %...finally add this subject's data into a row of PLS
                        %matrix:
                        PLSmat.(Measures{iM}){iG}(RowIND,:) = temp(:);
                        
                    end  %for Nm
                    clear MJLD;
                    %...remove the used file...
                    SubjsFiles{1,iG}{RowIND,1}=D(iF).name;
                    D(iF)=[];
                    Nf=numel(D)-2;
                    break
                    
                end% if choose file
                
                 catch
                        keyboard;
                end
                
            end %for Nf
        end%for Ns
        
        %...calculate also statistics (mean and std) of all subjects for
        %this condition and group..
        for iM = 1:Nm;%...for each measure...
            RowIND = (iC-1)*Ns(iG)+[1:Ns(iG)];
            PLSmatStats.(Measures{iM}).mean{iG,iC} = mean( PLSmat.(Measures{iM}){iG}(RowIND,:) );
            PLSmatStats.(Measures{iM}).std{iG,iC} = std( PLSmat.(Measures{iM}){iG}(RowIND,:) );
            
            %...reshape to original sizes of the measures...
            PLSmatStats.(Measures{iM}).mean{iG,iC} = reshape( PLSmatStats.(Measures{iM}).mean{iG,iC}, MeasOutsStruct{iM} );
            PLSmatStats.(Measures{iM}).std{iG,iC} = reshape( PLSmatStats.(Measures{iM}).std{iG,iC}, MeasOutsStruct{iM} );
            
            temp = PLSmat.(Measures{iM}){iG}(RowIND,:);
            ind = find(isnan(temp));
            if ~isempty(ind)
                warning('NaN values in PLS matrix converted to 0s')
                iG
                Groups{iG}
                iC
                Conds{iC}
                iS
                iF
                D(iF).name
                iM
                Measures{iM}
                temp(ind) = 0;
                PLSmat.(Measures{iM}){iG}(RowIND,:) = temp;
            end
            
        end
    end%for Nc
end %for Ng

cfgMeasure=cfgPS;
cfgMeasure.freqs = cfgPS.f;
cfgMeasure.scales = cfgPS.NchPairs;
save([PLSmatFolder,filesep,'PLSmat.mat'],'PLSmat','PLSmatStats','Groups','Conds','Subjs','SubjsFiles','Measures','MeasOutsStruct',...
    'Ng','Nc','Ns','Nm','Ntr')
save([StatsFolder,filesep,'Stats.mat'],'PLSmatStats','Groups','Conds','Subjs','SubjsFiles','Measures','MeasOutsStruct',...
    'Ng','Nc','Ns','Nm','Ntr')
if exist('cfg','var')
    fields = fieldnames(cfg);
    for i=1:numel(fields);
        cfgMeasure.(fields{i}) = cfg.(fields{i});
    end
end
if exist('cfgMeasure','var')
    save([PLSmatFolder,filesep,'PLSmat.mat'],'cfgMeasure','-append')
end
if exist('Colors','var')
    save([PLSmatFolder,filesep,'PLSmat.mat'],'Colors','-append')
end

%---------------------For Plotting-----------------------------------------
% for iM=1:Nm;
%     h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,AllMeasures, Colors, Groups, Conds, Ns,cfgMeasure);
%     for iiM = 1:length(AllMeasures);
%         fileName = [StatsFolder, filesep, AllMeasures{iM},'_statsMeanSTDerr.fig'];
%         saveas(h(iiM),fileName);
%     end
%     close(h);
%     clear h;
% end

for iM=1:Nm;
    h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM}, Colors, Groups, Conds, [],cfgMeasure);%AllMeasures
%     for iiM = 1:length(AllMeasures);
        fileName = [StatsFolder, filesep,Measures{iM},'_statsMean.fig'];% AllMeasures{iiM}
        saveas(h,fileName); %(iiM)
%     end
    close(h);
    clear h;
end
%---------------------For Plotting-----------------------------------------

toc

%clear all;