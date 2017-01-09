%function [PLSres, PLScfg] = DPrunPLStask( PLSfolder, RESfolder,  RESstruct, Conds, Groups, Measures, MeasOutsStruct, Ntr, Ns, Subjs)
%DPrunPLS This function calculates PLS analysis for several measures on a
%potentially multigroup and multisubject dataset

%   Outputs:
% -PLSres: structure with PLS results for each measure
% -PLScfg: structure with PLS configuration (PLSmat included)  for each measure

%% %---------------------A. Stack PLS matrix----------------------------------
%   Inputs:
RESfolder='./Results';
%-RESfolder: full path to the folder of results .mat files
PLSmatFolder = 'C:\Users\dionperd\Dropbox\Dionysis\PostDocMRS\Results\EEG_PSD_JLD\PLSallParamsTogether';
StatsFolder = [PLSmatFolder,filesep,'Stats'];
if ~exist(StatsFolder,'dir')
    mkdir(StatsFolder)
end
%-PLSmatFolder: full path to the folder where PLS matrix will be saved
RESstruct = 'MJLD';
%-RESstruct: name of the result structure (string)
Conds={'RestEC','Odd_nC','Odd_C'};%'RestEC','Odd_nC','Odd_C'
%-Conds: cellstring of the names of conditions
Groups={'mdo11','mdo12','mdo13','mdo14'};% 
%-Groups: cellstring of the names of Groups 
%(a string common to all files, e.g., the extension, if there is only 1 group)
Ns= [24 28 31 28]; %   % RestEO
%Ns= [30 26]; %  23 28 % The rest conditions without outliers
%-Ns: number of subjects of each group (vector of positive integers)
%load('C:\Users\CoordAgeEEG\Dropbox\Results\Subjects.mat')
Subjs = {};%
%-Subj: cell of cellstrings of the names of subjects for each group in case we want to use a
%       specific subset of them, otherwise Subjs={}
Measures = {'paramsA'};%{'alpha','beta','sigma','mu','RMSEcdf'};
%-Measures: a cellstring of the names of measures to be calculated
%                   all 4 JLD params   
%                  scales>1.004sec 10 freqs 4 params  
%MeasOutsStruct = {[20              10      4 ]};%{[76,10],[76,10],[76,10],[76,10],[76, 10]};
MeasOutsStruct =  {[76              10      4 ]};%{[76,10],[76,10],[76,10],[76,10],[76, 10]};

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
PlotGroups={'YC','OC','YA','OA'};%
PlotConds={'REC','UOT','AOT'}; %'RestEO',
%Colors = {rgb('DodgerBlue'),rgb('ForestGreen'),rgb('Red'), rgb('SteelBlue'),rgb('Teal'),rgb('DarkRed')};
Colors = {...
          rgb('Magenta'),rgb('DarkOrchid'),rgb('Indigo'),...
          rgb('Crimson'),rgb('IndianRed'),rgb('Maroon')...
          rgb('Lime'),rgb('ForestGreen'),rgb('DarkGreen'),...
          rgb('Cyan'),rgb('SteelBlue'),rgb('DarkBlue')...
          };
AllMeasures = {'\alpha','\beta','\sigma','\mu'};
%varargin={}; %any other inputs you would like to give to the plotting function
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds, Ns)plotMeanSTD_2D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
%plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns)plotMeanSTD_1D(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,varargin);
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds)plotMeanSTD_0D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
%method = 'MeanStdErr';
%%plotMeasure = @(M,STD,Ns,Ng,Nc,Groups,Conds,method,measure,cfg)plotUnivComplmapMeanSTD(M,STD,Ns,Ng,Nc,PlotGroups,PlotConds,method,measure,cfg);
%plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns,scales)plotMeanSTD_1D_JLDparams(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,scales);
plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns,scales,freqs)plotMeanSTD_2D_JLDparamsAlltogether(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,cfgMeasure);
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
Nf=numel(D)-2;
%The function to choose among the files...
if (isempty(Subjs))
    %...without subjects selection...
    ChooseFileCommand = '( ~isempty(strfind(D(iF+2).name,Groups{iG})) && ~isempty(strfind(D(iF+2).name,Conds{iC})) )';
else
    %...and with subjects selection...
    ChooseFileCommand = '( ~isempty(strfind(D(iF+2).name,Groups{iG})) && ~isempty(strfind(D(iF+2).name,Conds{iC}))&& ~isempty(strfind(D(iF+2).name,Subjs{iG}{iS})) )';
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
                    load(D(iF+2).name)
                    
                    disp('PLSmat row index: ')
                    RowIND = (iC-1)*Ns(iG)+iS;
                    disp(RowIND)
                    
                    JLD.alpha = squeeze(JLD.paramsA(1,:,:));
                    JLD.alpha = DPcheckMeasureForPLSmat(JLD.alpha,0,2,D(iF+2).name,Groups{iG},Conds{iC},'alpha',iG,iC,iS,iF);

                    JLD.beta = squeeze(JLD.paramsA(2,:,:));
                    JLD.beta = DPcheckMeasureForPLSmat(JLD.beta,-1,1,D(iF+2).name,Groups{iG},Conds{iC},'beta',iG,iC,iS,iF);
                    
                    JLD.sigma = squeeze(JLD.paramsA(3,:,:));
                    JLD.sigma = DPcheckMeasureForPLSmat(JLD.sigma,-1,1,D(iF+2).name,Groups{iG},Conds{iC},'sigma',iG,iC,iS,iF);
                    
                    JLD.mu = squeeze(JLD.paramsA(4,:,:));
                    JLD.mu = DPcheckMeasureForPLSmat(JLD.mu,-1,1,D(iF+2).name,Groups{iG},Conds{iC},'mu',iG,iC,iS,iF);
                    
%                    JLD.RMSEcdf = JLD.RMSElogcdf;
                    JLD = rmfield(JLD,'paramsA');
                    JLD.paramsA(:,:,1) = JLD.alpha;
                    JLD.paramsA(:,:,2) = JLD.beta;
                    JLD.paramsA(:,:,3) = JLD.sigma;
                    JLD.paramsA(:,:,4) = JLD.mu;
                     
                    for iM = 1:Nm;%...for each measure...
                        
                        %...get this subject's data into temp...
                        temp=eval(tempCommand{iM});
                        
                        %                         disp(Measures{iM})
                        %                         disp(size(temp))
                        
                        %...finally add this subject's data into a row of PLS
                        %matrix:
                        PLSmat.(Measures{iM}){iG}(RowIND,:) = temp(:);
                        
                    end  %for Nm
                    clear JLD;
                    %...remove the used file...
                    SubjsFiles{1,iG}{RowIND,1}=D(iF+2).name;
                    D(iF+2)=[];
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
                Conds{iG}
                iS
                iF
                D(iF+2).name
                iM
                Measures{iM}
                temp(ind) = 0;
                PLSmat.(Measures{iM}){iG}(RowIND,:) = temp;
            end
            
        end
    end%for Nc
end %for Ng

cfgMeasure.scales = cfgJLD.scales*cfgPS.Ts;
cfgMeasure.freqs = cfgJLD.f;
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
%     h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,AllMeasures, Colors, Groups, Conds, Ns);
%     for iiM = 1:length(AllMeasures);
%         fileName = [StatsFolder, filesep, AllMeasures{iM},'_statsMeanSTDerr.fig'];
%         saveas(h(iiM),fileName);
%     end
%     close(h);
%     clear h;
% end

for iM=1:Nm;
    h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,AllMeasures, Colors, Groups, Conds, []);
    for iiM = 1:length(AllMeasures);
        fileName = [StatsFolder, filesep, AllMeasures{iiM},'_statsMean.fig'];
        saveas(h(iiM),fileName);
    end
    close(h);
    clear h;
end
%---------------------For Plotting-----------------------------------------

toc

%clear all;