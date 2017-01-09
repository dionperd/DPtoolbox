%function [PLSres, PLScfg] = DPrunPLStask( PLSfolder, RESfolder,  RESstruct, Conds, Groups, Measures, MeasOutsStruct, Ntr, Ns, Subjs)
%DPrunPLS This function calculates PLS analysis for several measures on a
%potentially multigroup and multisubject dataset

%   Outputs:
% -PLSres: structure with PLS results for each measure
% -PLScfg: structure with PLS configuration (PLSmat included)  for each measure

%% %---------------------A. Stack PLS matrix----------------------------------
%   Inputs:
RESfolder='C:\Users\dionperd\Dropbox\Dionysis\PostDocMRS\Results\EEG_PSD_JLD\ResultsFiles';
%-RESfolder: full path to the folder of results .mat files
PLSmatFolder = 'C:\Users\dionperd\Dropbox\Dionysis\PostDocMRS\Results\EEG_PSD_JLD\PLS\';
StatsFolder = [PLSmatFolder,filesep,'Stats'];
if ~exist(StatsFolder,'dir')
    mkdir(StatsFolder)
end
%-PLSmatFolder: full path to the folder where PLS matrix will be saved
RESstruct = 'JLD';
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
Measures = {'alpha','beta','sigma','mu','RMSEcdf'};%{'P','PS','DOF','logF','H','logV','STD','MSE','MSEn'};%'P','DOF','PS','logV','STD','MSE','logF','H'
%-Measures: a cellstring of the names of measures to be calculated
%                                                                                       P       PS    DOF    logF     H     logV     
MeasOutsStruct = {[76,10],[76,10],[76,10],[76,10],[76, 10]};%{[2048 58],[1 58],[1 58],[47 58],[1 58],[50 58],...
                                                                                        %[50 58],[50 58],[50 58]};
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
PlotConds={'R','OnC','OC'}; %'RestEO',
%Colors = {rgb('DodgerBlue'),rgb('ForestGreen'),rgb('Red'), rgb('SteelBlue'),rgb('Teal'),rgb('DarkRed')};
Colors = {...
          rgb('Magenta'),rgb('DarkOrchid'),rgb('Indigo'),...
          rgb('Crimson'),rgb('IndianRed'),rgb('Maroon')...
          rgb('Lime'),rgb('ForestGreen'),rgb('DarkGreen'),...
          rgb('Cyan'),rgb('SteelBlue'),rgb('DarkBlue')...
          };
%varargin={}; %any other inputs you would like to give to the plotting function
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds, Ns)plotMeanSTD_2D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
%plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns)plotMeanSTD_1D(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,varargin);
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds)plotMeanSTD_0D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
%method = 'MeanStdErr';
%%plotMeasure = @(M,STD,Ns,Ng,Nc,Groups,Conds,method,measure,cfg)plotUnivComplmapMeanSTD(M,STD,Ns,Ng,Nc,PlotGroups,PlotConds,method,measure,cfg);
%plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns,scales)plotMeanSTD_1D_JLDparams(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,scales);
scales = 0.004*unique(round(exp(log(1):(log(1000)-log(1))/100:log(1000))));
freqs = [2:2:20];
cfgMeasure.scales = scales;
cfgMeasure.freqs = freqs;
plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns,scales,freqs)plotMeanSTD_2D_JLDparams(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,cfgMeasure);
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
                
                
                %...if this file corresponds to this group, condition (and subject)...
                if eval(ChooseFileCommand)

                    %...load the file...
                    load(D(iF+2).name)
                    
                    disp('PLSmat row index: ')
                    RowIND = (iC-1)*Ns(iG)+iS;
                    disp(RowIND)
                    
                    JLD.alpha = squeeze(JLD.paramsA(1,:,:));
                    JLD.beta = squeeze(JLD.paramsA(2,:,:));
                    JLD.sigma = squeeze(JLD.paramsA(3,:,:));
                    JLD.mu = squeeze(JLD.paramsA(4,:,:));
                    JLD.RMSEcdf = JLD.RMSElogcdf;
                    
                    for iM = 1:Nm;%...for each measure...
                        
                        %...get this subject's data into temp...
                        temp=eval(tempCommand{iM});
                        
                        %                         disp(Measures{iM})
                        %                         disp(size(temp))
                        
                        %...finally add this subject's data into a row of PLS
                        %matrix:
                        PLSmat.(Measures{iM}){iG}(RowIND,:) = temp(:);
                        
                    end  %for Nm
                    
                    %...remove the used file...
                    SubjsFiles{1,iG}{RowIND,1}=D(iF+2).name;
                    D(iF+2)=[];
                    Nf=numel(D)-2;
                    break
                end% if choose file
            end %for Nf
        end%for Ns
        
        %...calculate also statistics (mean and std) of all subjects for
        %this condition and group..
        for iM = 1:Nm;%...for each measure...
            RowIND = (iC-1)*Ns(iG)+[1:Ns(iG)];
            PLSmatStats.(Measures{iM}).mean{iG,iC} = nanmean( PLSmat.(Measures{iM}){iG}(RowIND,:) );
            PLSmatStats.(Measures{iM}).std{iG,iC} = nanstd( PLSmat.(Measures{iM}){iG}(RowIND,:) );
            
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
for iM=1:Nm;
    h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM}, Colors, Groups, Conds, Ns,cfgMeasure);
    fileName = [StatsFolder, filesep, Measures{iM},'_statsMeanStdErr.fig'];
    saveas(h,fileName);
    close(h);
    clear h;
end
%---------------------For Plotting-----------------------------------------

%---------------------For Plotting-----------------------------------------
for iM=1:Nm;
    h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM}, Colors, Groups, Conds, [],cfgMeasure);
    fileName = [StatsFolder, filesep, Measures{iM},'_stats.fig'];
    saveas(h,fileName);
    close(h);
    clear h;
end
%---------------------For Plotting-----------------------------------------

toc

clear all;





%% %-------------------------B. Run PLS--------------------------------------
%Inputs:
PLSmatFolder = 'C:\Users\dionperd\Dropbox\Dionysis\PostDocMRS\Results\EEG_PSD_JLD\PLS';
%-PLSmatFolder: full path to the folder where PLS matrix is saved
OutputFolder= 'C:\Users\dionperd\Dropbox\Dionysis\PostDocMRS\Results\EEG_PSD_JLD\PLS\MeanCntr';
if ~exist(OutputFolder,'dir')
    mkdir(OutputFolder)
end
%-OutputFolder: full path to the folder where PLS output will be stored
PLSfolder='C:\Users\dionperd\Dropbox\Dionysis\DPtoolbox\Statistics\external\plscmd';
%-PLSfolder: full path to the plscmd folder

%-Options structure opt for PLS
opt.method = 1;
%	method: This option will decide which PLS method that the
%		program will use:
%			1. Mean-Centering Task PLS
%			2. Non-Rotated Task PLS
%			3. Regular Behavior PLS
%			5. Non-Rotated Behavior PLS
%		If it is not specified, program will use default value 1.
%
%--------------------------------------------------------------------------
opt.num_perm = 1000;
%	num_perm: If specified and num_perm > 0, PLS will run permutation
%		test with num_perm amount of samples; otherwise, num_perm
%		will use default value 0, and program will not run 
%		permutation test.
%
%	is_struct: If it is not specified, is_struct will use default
%		value 0. If specified and is_struct = 1, PLS will not 
%		permute conditions within a group. You only need to 
%		specify it when you run Non-Behavior Structure PLS and
%		the segmented data is used.
%
%	num_split: If specified and num_split > 0, PLS will run permutation
%		test with Natasha's Split Half routine; otherwise num_split
%		will use default value 0, and program will not use Split 
%		Half routine.
%
%--------------------------------------------------------------------------
opt.num_boot = 1000;
%	num_boot: If specified and num_boot > 0, PLS will run bootstrap
%		test with num_boot amount of samples; otherwise, num_boot
%		will use default value 0, and program will not run
%		bootstrap test.
%
%--------------------------------------------------------------------------
%opt.clim = 99;
%	clim: Confidence level between 0 and 100. If not specified,
%		program will use default value 95.
%
%--------------------------------------------------------------------------
%load('path to designdata matrix')
%opt.stacked_designdata = DesignData;
%	stacked_designdata: If you are choosing Non-Rotated Task PLS,
%		you have to specify this 2-D numerical matrix.
%
%		Number of columns always stand for number of designs. 
%
%		For Non-Rotated Task PLS, number of rows in each group
%		equal to number of conditions. So total number
%		of rows equal to that multiplied by number of groups, in
%		the form of "condition in group". i.e.:
%			g1	c1
%				c2
%				c3
%			g2	c1
%				c2
%				c3
%
%		For Non-Rotated Behavior PLS, number of rows in each
%		group equal to number of conditions multiplied by number
%		of behavior measures. So total number of rows equal to
%		that multiplied by number of groups, in the form of
%		"measure in condition in group". i.e.:
%			g1	c1	b1
%					b2
%				c2	b1
%					b2
%				c3	b1
%					b2
%			g2	c1	b1
%					b2
%				c2	b1
%					b2
%				c3	b1
%					b2
%
Nb = 0; %initialize number of behavioral measures
% load('C:\data\DATA_Aging_Complexity\ResultsFinal\');
% opt.stacked_behavdata = Age;
% Nb = size(BehavData,2); %number of behavioral measures
% Behav = {'ReactionTime'}; %cell of strings with behavioral measures'
% labels
%	stacked_behavdata: If you are choosing any Behavior PLS or
%		Multiblock PLS, you have to specify this 2-D numerical
%		matrix.
%
%		Number of columns always stand for behavior measures.
%		Number of rows should equal to the sum of the number
%		of rows in datamat list array, and it is in the form
%		of "subject in condition in group".
%
%		This option can be applied to method 3, 4, 5 and 6.
%--------------------------------------------------------------------------
opt.meancentering_type=2;
%	meancentering_type: Type of meancentering.
%		0. Remove group condition means from conditon means
%		   within each group. Tells us how condition effects
%		   are modulated by group membership. (Boost condition
%		   differences, remove overall group differences).
%		1. Remove grand condition means from each group condition
%		   mean. Tells us how conditions are modulated by group
%		   membership (Boost group differences, remove overall
%		   condition differences).
%		2. Remove grand mean over all subjects and conditions.
%		   Tells us full spectrum of condition and group effects.
%		3. Remove all main effects by subtracting condition and
%		   group means. This type of analysis will deal with
%		   pure group by condition interaction.
%		If it is not specified, program will use default value 0.
%		This option can be applied to method 1, 2, 4, and 6.
%
%--------------------------------------------------------------------------
%	cormode: correaltion mode determines correlation type to analyze.
%		0. Pearson correlation
%		2. covariance
%		4. cosine angle
%		6. dot product
%
%		If it is not specified, program will use default value 0.
%		This option can be applied to method 3, 4, 5, and 6.
%--------------------------------------------------------------------------
%	boot_type: 'strat' (default for standard PLS approach), or
%		'nonstrat' for nonstratified boot samples.
%----------------------------End of Inputs---------------------------------

cd(PLSmatFolder) 
load('PLSmat.mat')
% cfgMeasure.scales = cfgMeasure.scales(cfgMeasure.scales>2.08); 
% ind = [];
% for ii=1:10;
%     ind = [ind,76*(ii-1)+[67:76]];
% end
%Run PLS
addpath(PLSfolder); %add PLS cmd path
disp('Running PLS for measure:...')
for iM = 1:Nm;%...for each measure...
    
    disp(['...',Measures{iM},'...'])
    tic
    %Run PLS
    %result = pls_analysis(datamat_lst, num_subj_lst, num_cond, option);
%     for iG = 1:Ng;
%         PLSmat.(Measures{iM}){iG} = PLSmat.(Measures{iM}){iG}(:,ind); 
%     end
   PLSres.(Measures{iM}) = DPpls_analysis(PLSmat.(Measures{iM}), Ns, Nc, opt);
   toc
    
    %Store PLSmat and Msiz(iM)
    PLScfg.(Measures{iM}).PLSmat = PLSmat.(Measures{iM});
    PLSmat = rmfield(PLSmat,Measures{iM});
    PLScfg.(Measures{iM}).MeasOutsStruct = MeasOutsStruct{iM}; %[10 10]
    
    
    %Create meaningful ouputs out of this result:
    
% Now this output is returned by DPpls_analysis by default for mean
% centered PLS
% %     %A. Percentage of variance explained by each Latent Variable
%     s2 = (PLSres.(Measures{iM}).s).^2; %square of eigenvalues
%     PLSres.(Measures{iM}).var100 = 100*s2/sum(s2);     
            
    
   %B. Brain scores
    Nlv = size(PLSres.(Measures{iM}).s,1);% the number of LVs   
    Ndof = size(PLSres.(Measures{iM}).v,1);%number of degrees of freedom plus 1 (groups x conditions x behavioral measures)
    Nr = size(PLSres.(Measures{iM}).usc,1); % the number of brain scores
    %reshape and normalize brain scores
    %uscNorm = PLSres.(Measures{iM}).usc./repmat((PLSres.(Measures{iM}).s).',[Nr,1]);
    %separate brain scores for each group and condition

    for iLV = 1:Nlv;
        for iG = 1:Ng;
            for iC = 1:Nc;
                PLSres.(Measures{iM}).BrainScores{iLV}{iG,iC} = PLSres.(Measures{iM}).usc(Nc*sum(Ns(1:iG-1))+(iC-1)*Ns(iG)+[1:Ns(iG)],iLV);
                PLSres.(Measures{iM}).BrainScoresMU{iLV}(iG,iC) = mean(PLSres.(Measures{iM}).BrainScores{iLV}{iG,iC});
                PLSres.(Measures{iM}).BrainScoresSTDERR{iLV}(iG,iC) = std(PLSres.(Measures{iM}).BrainScores{iLV}{iG,iC})/sqrt(Ns(iG));
            end
        end
        PLSres.(Measures{iM}).BrainScoresGrandMean(iLV) = mean(PLSres.(Measures{iM}).usc(:,iLV));
    end
    
    %B. Brain LVs 
    for iLV=1:Nlv;
            PLSres.(Measures{iM}).BrainLVs{iLV} = reshape(PLSres.(Measures{iM}).u(:,iLV),PLScfg.(Measures{iM}).MeasOutsStruct);
            PLSres.(Measures{iM}).BrainLVsBootsRatios{iLV} = [];
    end 
    
    %...if there is a bootstrap result (get also brain and task LVs' standard errors)...
    if isfield(PLSres.(Measures{iM}),'boot_result')
 
        %reshape the boostrap standard error and error corrected brain scores
        %PLSres.(Measures{iM}).boot_result.compare_u  = reshape(PLSres.(Measures{iM}).boot_result.compare_u,[PLSres.(Measures{iM}).MeasOutsStruct, Nlv]);
        %PLSres.(Measures{iM}).boot_result.u_se  = reshape(PLSres.(Measures{iM}).boot_result.u_se,[PLScfg.(Measures{iM}).MeasOutsStruct, Nlv]);
        
        for iLV = 1:Nlv;
                    
            PLSres.(Measures{iM}).BrainLVsBootsRatios{iLV} = reshape(PLSres.(Measures{iM}).boot_result.compare_u(:,iLV),[PLScfg.(Measures{iM}).MeasOutsStruct]);

            for iG = 1:Ng;
                for iC = 1:Nc;
                    PLSres.(Measures{iM}).BrainScoresMeanCntr{iLV}{iG,iC} = PLSres.(Measures{iM}).boot_result.usc2(Nc*sum(Ns(1:iG-1))+(iC-1)*Ns(iG)+[1:Ns(iG)],iLV);
                end
            end
            PLSres.(Measures{iM}).BrainScoresBootsMeanCntrMU{iLV} = reshape( PLSres.(Measures{iM}).boot_result.orig_usc(:,iLV) ,Ng,Nc );
            %if any(isnan(PLSres.(Measures{iM}).boot_result.llusc_adj)) || any(isnan(PLSres.(Measures{iM}).boot_result.ulusc_adj))
                PLSres.(Measures{iM}).BrainScoresBootsLL{iLV} =  reshape( PLSres.(Measures{iM}).boot_result.llusc(:,iLV) ,Ng,Nc );
                PLSres.(Measures{iM}).BrainScoresBootsUL{iLV} =  reshape( PLSres.(Measures{iM}).boot_result.ulusc(:,iLV) ,Ng,Nc );
            %else
            %    PLSres.(Measures{iM}).BrainScoresBootsLL{iLV} =  reshape( PLSres.(Measures{iM}).boot_result.llusc_adj(:,iLV) ,Ng,Nc );
            %    PLSres.(Measures{iM}).BrainScoresBootsUL{iLV} =  reshape( PLSres.(Measures{iM}).boot_result.ulusc_adj(:,iLV) ,Ng,Nc );
           % end
        end

    end  
end
PLScfg.opt = opt;
PLScfg.Ng = Ng;
PLScfg.Nc = Nc;
PLScfg.Nb = Nb;
PLScfg.Ns = Ns;
PLScfg.Nm = Nm;
PLScfg.Nr = Nr;
PLScfg.Ntr = Ntr;
PLScfg.Nlv =Nlv;
PLScfg.Ndof = Ndof;
PLScfg.Conds = Conds;
PLScfg.Groups = Groups;
if opt.method==3
    PLScfg.Behav = Behav;
end
PLScfg.Subjs = Subjs;
PLScfg.SubjsFiles = SubjsFiles;
PLScfg.Measures=Measures;
PLScfg.MeasOutsStruct=MeasOutsStruct;
if exist('cfgMeasure','var')
    PLScfg.cfgMeasure = cfgMeasure;
end
if exist('Colors','var')
    PLScfg.Colors = Colors;
end
disp('Saving results...')
filePathName = fullfile(OutputFolder,'PLSres.mat');
save(filePathName,'PLSres','PLScfg','OutputFolder')
disp('DONE!')
clear all;





%% %-------------------------C. Plotting PLS---------------------------------
%Inputs:
OutputFolder= 'C:\Users\dionperd\Dropbox\Dionysis\PostDocMRS\Results\EEG_PSD_JLD\PLS\MeanCntr';
%-OutputFolder: full path to the folder where PLS output will be stored
alpha = 0.1; %statistical alpha value for significance
%-alpha: alpha value for significance of permutation test of LVs (real 0<scalar<1)
%varMin = 1; %minimum % variance explained by a LV to be plotted
%Method for plotting brain LVs:
%-'BootsRatios': plot bootstrap ratios (possible only if there is a
%                bootstrap test done)
%-'LV': plot the Brain LV
BrainLVplotOpts.method = 'BootsRatios';
%threshold for reliability
BrainLVplotOpts.th = 0; %~99% confidence
%alpha range for points below threshold for maps of 2D Brain LVs
BrainLVplotOpts.Range = [0 0.5]; 
%-alpha: alpha value for significance of permutation test of LVs (real 0<scalar<1)
%Behav = {'Age'};

cd(OutputFolder) 
load('PLSres.mat')
Ng = PLScfg.Ng;
Nc = PLScfg.Nc;

Conds={'R','OnC','OC'};%'RestEC','Odd_nC','Odd_C'
%-Conds: cellstring of the names of conditions
Groups={'YC','OC','YA','OA'};% 
%-Groups: cellstring of the names of Groups

if isfield(PLScfg,'cfgMeasure')
    cfgMeasure = PLScfg.cfgMeasure;
end

if isfield(PLScfg,'Colors')
    Colors = PLScfg.Colors;
else
    Colors = {rgb('DodgerBlue'),rgb('ForestGreen'),rgb('Red'), rgb('DarkBlue'),rgb('Teal'),rgb('DarkRed')};
end
%plotBrainLVs = @(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,cfgMeasure)plotUnivComplmap(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,cfgMeasure); 
plotBrainLVs = @(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,scales,freqs)plotBrainLVs_2D_JLDparams(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,cfgMeasure);
% plotBrainLVs = @(BrainLV,Measure,cfgMeasure)plotBrainLVs_1D(BrainLV,BrainLVbootsRatio,BrainLVplotOpts,Measure,cfgMeasure,varargin);
%-plotBrainLVs: a handle to a function to plot the brain LVs of each measure
%----------------------------End of Inputs---------------------------------


iX=0;
if PLScfg.opt.method<3
    for iG=1:PLScfg.Ng;
        for iC=1:PLScfg.Nc;
            iX=iX+1;
            if Ng>1
                xlbl{iX} = [Groups{iG},'-',Conds{iC}];
            else
                xlbl{iX} = Conds{iC};
            end
        end
    end
else
    Nb = PLScfg.Nb;
    for iG=1:PLScfg.Ng;
        for iC=1:PLScfg.Nc;
            for iB=1:PLScfg.Nb;
                iX=iX+1;
                if Ng>1
                    xlbl{iX} = [Groups{iG},'-',Conds{iC},'-',Behav{iB}];
                else
                    xlbl{iX} = [Conds{iC},'-',Behav{iB}];
                end    
            end
        end
    end
end
Measures = PLScfg.Measures;
Nlv = PLScfg.Nlv;
Ndof = PLScfg.Ndof;
disp('Plotting PLS results for measure:...')

for iM = 1:PLScfg.Nm;%...for each measure...
    
    thisPLS = PLSres.(PLScfg.Measures{iM});

    NlvSign = Nlv;
    LVsign = 1:NlvSign;
    if isfield(thisPLS,'perm_result') 
        %Choose the LVs that have p<= alpha and explain more than varMin %
        %of variance:
        indLVsign = (thisPLS.perm_result.sprob <= alpha);% & (thisPLS.var100 >= varMin);
        LVsign = LVsign(indLVsign);
        NlvSign = length(LVsign);
    end
    
    
    h=figure('Name',[PLScfg.Measures{iM}, ', task saliences and brain scores'],'units','normalized','outerposition',[0 0 1 1]);

    for iLV = 1:NlvSign;
        
        LV  = LVsign(iLV);
        
        
        if thisPLS.v(1,iLV) > 0
            sgn = -1;
        else
            sgn = 1;
        end
        
        figure(h);
        
        subplot(3,NlvSign,iLV)
        if any(PLScfg.opt.method==[1,3]) %if the method is mean-centered...
            title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, ',num2str(thisPLS.var100(LV)),'%, p = ',num2str(thisPLS.perm_result.sprob(LV)),'):'];['task LVs and standardized brain scores'' distribution']})
        else %if it is non-rotated...
            title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, s = ',num2str(thisPLS.s(LV)),', p = ',num2str(thisPLS.perm_result.sprob(LV)),'):'];['task LVs and standardized brain scores'' distribution']})
        end
        hold on;
        for iX=1:Ndof;
            bar(iX,sgn*thisPLS.v(iX,LV),0.5,'FaceColor',Colors{iX});
            hold on;
            if isfield(PLSres.(Measures{iM}),'boot_result')
                errorbar(iX,sgn*thisPLS.BrainScoresBootsMeanCntrMU{iLV}(iX)/thisPLS.s(iLV),...
                         min([sgn*thisPLS.BrainScoresBootsLL{iLV}(iX),sgn*thisPLS.BrainScoresBootsUL{iLV}(iX)])/thisPLS.s(iLV),...
                         max([sgn*thisPLS.BrainScoresBootsLL{iLV}(iX),sgn*thisPLS.BrainScoresBootsUL{iLV}(iX)])/thisPLS.s(iLV),'k')
                plot(iX,sgn*thisPLS.BrainScoresBootsMeanCntrMU{iLV}(iX)/thisPLS.s(iLV),'k*','Markersize',5)
            else
                errorbar(iX,sgn*thisPLS.BrainScoresMU{iLV}(iX)-sgn*thisPLS.BrainScoresGrandMean(iLV),...
                         (sgn*thisPLS.BrainScoresMU{iLV}(iX)-sgn*thisPLS.BrainScoresGrandMean(iLV)...
                         -thisPLS.BrainScoresSTDERR{iLV}(iX))/thisPLS.s(iLV),...
                         (sgn*thisPLS.(Measures{iM}).BrainScoresMU{iLV}(iX)-sgn*thisPLS.BrainScoresGrandMean(iLV)...
                         +thisPLS.BrainScoresSTDERR{iLV}(iX))/thisPLS.s(iLV),'k')
                plot(iX,sgn*thisPLS.BrainScoresMU{iLV}(iX)-sgn*thisPLS.BrainScoresGrandMean(iLV),'k*','Markersize',5)
            end
        end
        set(gca,'xlim',[0 Ndof+1],'xtick',[1:Ndof],'xticklabel',xlbl);
        for iX = 1:Ndof;
            if numel(xlbl{iX})>5
                xticklabel_rotate([],45)
                break;
            end
        end
        hold off
        
        subplot(3,NlvSign,NlvSign+iLV)
        if any(PLScfg.opt.method==[1,3]) %if the method is mean-centered...
            title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, ',num2str(thisPLS.var100(LV)),'%, p = ',num2str(thisPLS.perm_result.sprob(LV)),'):'];['brain scores']})
        else %if it is non-rotated...
            title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, s = ',num2str(thisPLS.s(LV)),', p = ',num2str(thisPLS.perm_result.sprob(LV)),'):'];['brain scores']})
        end
        for iX=1:Ndof;
            hold on;
            if isfield(PLSres.(Measures{iM}),'boot_result')
                thisM = (sgn*thisPLS.BrainScoresBootsMeanCntrMU{iLV}(iX)+sgn*thisPLS.BrainScoresGrandMean(iLV));
                bar(iX,thisM,0.5,'FaceColor',Colors{iX});
                errorbar(iX,thisM,...
                         min([sgn*thisPLS.BrainScoresBootsLL{iLV}(iX),sgn*thisPLS.BrainScoresBootsUL{iLV}(iX)]),...
                         max([sgn*thisPLS.BrainScoresBootsLL{iLV}(iX),sgn*thisPLS.BrainScoresBootsUL{iLV}(iX)]),'k')
            else
                bar(iX,sgn*thisPLS.BrainScoresMU{iLV}(iX),0.5,'FaceColor',Colors{iX});
                errorbar(iX,sgn*thisPLS.BrainScoresMU{iLV}(iX),...
                         sgn*thisPLS.BrainScoresMU{iLV}(iX) - thisPLS.BrainScoresSTDERR{iLV}(iX),...
                         sgn*thisPLS.(Measures{iM}).BrainScoresMU{iLV}(iX) + thisPLS.BrainScoresSTDERR{iLV}(iX),'k')
            end
        end
        set(gca,'xlim',[0 Ndof+1],'xtick',[1:Ndof],'xticklabel',xlbl);
        for iX = 1:Ndof;
            if numel(xlbl{iX})>5
                xticklabel_rotate([],45)
                break;
            end
        end
        hold off
        
        subplot(3,NlvSign,2*NlvSign+iLV)
        if any(PLScfg.opt.method==[1,3]) %if the method is mean-centered...
            title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, ',num2str(thisPLS.var100(LV)),'%, p = ',num2str(thisPLS.perm_result.sprob(LV)),'):'];['brain scores per subject']})
        else %if it is non-rotated...
            title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, s = ',num2str(thisPLS.s(LV)),', p = ',num2str(thisPLS.perm_result.sprob(LV)),'):'];['brain scores per subject']})
        end
        hold on;
        iX=0;
        minVal = +inf;
        maxVal =-inf;
        for iG=1:PLScfg.Ng;
            for iC=1:PLScfg.Nc;
                iX=iX+1;
                if isfield(PLSres.(Measures{iM}),'boot_result') 
                    temp = double(sgn*thisPLS.BrainScoresMeanCntr{iLV}{iG,iC} + sgn*thisPLS.BrainScoresGrandMean(iLV));
                else
                    temp = double(sgn*thisPLS.BrainScores{iLV}{iG,iC});
                end
                minTemp=min(temp);
                if minTemp<minVal
                    minVal=minTemp;
                end
                maxTemp=max(temp);
                if maxTemp>maxVal
                    maxVal=maxTemp;
                end
                for iS = 1:PLScfg.Ns(iG);
                    text(iX,temp(iS),num2str(iS),'Color',Colors{iX});
                end
            end
        end
        set(gca,'xlim',[0 Ndof+1],'ylim',[minVal maxVal],'xtick',[1:Ndof],'xticklabel',xlbl);
        for iX = 1:Ndof;
            if numel(xlbl{iX})>5
                xticklabel_rotate([],45)
                break;
            end
        end
        hold off
         
        hb=figure('Name',[PLScfg.Measures{iM},', LV',num2str(LV),': brain latent variable'],'units','normalized','outerposition',[0 0 1 1]);
        if any(PLScfg.opt.method==[1,3]) %if the method is mean-centered...
            title([PLScfg.Measures{iM},', brain LV_{',num2str(LV),'} (',num2str(thisPLS.var100(LV)),'%, p = ',num2str(thisPLS.perm_result.sprob(LV)),')'])
        else
            title([PLScfg.Measures{iM},', brain LV_{',num2str(LV),'} (s = ',num2str(thisPLS.s(LV)),', p = ',num2str(thisPLS.perm_result.sprob(LV)),')'])
        end
        hold on;
        plotBrainLVs(sgn*thisPLS.BrainLVs{LV},sgn*thisPLS.BrainLVsBootsRatios{LV},BrainLVplotOpts,PLScfg.Measures{iM},cfgMeasure);
        hold off;
        fileNameB = [PLScfg.Measures{iM},'_LV',num2str(LV),'_BrainLV'];
        saveas(hb,[fileNameB,'.fig']);
        close(hb);
        
    end
    fileName = [PLScfg.Measures{iM},'_TaskLVs_BrainScores'];
    saveas(h,[fileName,'.fig']);
    close(h);
end