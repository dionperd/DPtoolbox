%function [PLSres, PLScfg] = DPrunPLStask( PLSfolder, RESfolder,  RESstruct, Conds, Groups, Measures, MeasOutsStruct, Ntr, Ns, Subjs)
%DPrunPLS This function calculates PLS analysis for several measures on a
%potentially multigroup and multisubject dataset

%Outputs:
%-PLSres: structure with PLS results for each measure
%-PLScfg: structure with PLS configuration (PLSmat included)  for each measure
scales(:,1) = 1;%unique(round(exp(log(1):(log(1000)-log(1))/100:log(1000))))/100;

% %---------------------A. Stack PLS matrix----------------------------------
%   Inputs:
RESfolder='C:\Users\vernooij.c\Documents\MATLAB\FCD\Results\PLS\JLD_paramsA';
%-RESfolder: full path to the folder of results .mat files
PLSmatFolder = 'C:\Users\vernooij.c\Documents\MATLAB\FCD\Results\PLS';
%-PLSmatFolder: full path to the folder where PLS matrix will be saved
RESstruct = 'JLD';
%-RESstruct: name of the result structure (string)
Conds={'ID1','ID2','ID3','ID4','ID5'};%{'RestEC','Odd_nC','Odd_C'};%'RestEC','Odd_nC','Odd_C'
%-Conds: cellstring of the names of conditions
Groups={'young','MSJLD_params'};% 'mdo11','mdo12',
%-Groups: cellstring of the names of Groups 
%(a string common to all files, e.g., the extension, if there is only 1 group)
Ns= [14 9]; % 24 28  % RestEO 
%-Ns: number of subjects of each group (vector of positive integers)
%load('C:\Users\CoordAgeEEG\Dropbox\Results\Subjects.mat')
Subjs = {};%
%-Subj: cell of cellstrings of the names of subjects for each group in case we want to use a
%       specific subset of them, otherwise Subjs={}
Measures = {'alpha','beta','sigma','mu'};%{'P','PS','DOF','logF','H','logV','STD','MSE','MSEn'};%'P','DOF','PS','logV','STD','MSE','logF','H'
%-Measures: a cellstring of the names of measures to be calculated
%                                                                                       P       PS    DOF    logF     H     logV     
MeasOutsStruct = {[1,1],[1,1],[1,1],[1,1]};%{[2048 58],[1 58],[1 58],[47 58],[1 58],[50 58],...
% MeasOutsStruct = {[4,1],[1,1],[1,1],[1,1]};                                                                            %[50 58],[50 58],[50 58]};
%-MeasOutsStruct: a cell of the same size as Measures with vectors
%                 corresponding to the size of the specific measure's 
%                 result PER TRIAL or for the mean data,
%                 but the trial index has to be the LAST one, if any
Ntr=0;
%-Ntr: number of trials to be considered for the within subject average (integer),
%      if Ntr = 0, an average over all trials is taken
%      if Ntr = 1, only the first trial will be taken, which is equivalent
%                  to taking mean data already calculated
%
%---------------------For Plotting-----------------------------------------
%Inputs for plotting:
PlotGroups={'Y','O'};%'YC','OC',
PlotConds={'ID1','ID2','ID3','ID4','ID5'};%;{'R','OnC','OC'}; %'RestEO',
%Colors = {rgb('DodgerBlue'),rgb('ForestGreen'),rgb('Red'), rgb('SteelBlue'),rgb('Teal'),rgb('DarkRed')};
Colors = {rgb('blue'),rgb('green'),rgb('yellow'),rgb('red'),rgb('black'),rgb('blue'),rgb('green'),rgb('yellow'),rgb('red'),rgb('black')};

% varargin={}; %any other inputs you would like to give to the plotting function
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds, Ns)plotMeanSTD_2D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
%plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns)plotMeanSTD_1D(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,varargin);
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds)plotMeanSTD_0D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
%method = 'MeanStdErr';
%%plotMeasure = @(M,STD,Ns,Ng,Nc,Groups,Conds,method,measure,cfg)plotUnivComplmapMeanSTD(M,STD,Ns,Ng,Nc,PlotGroups,PlotConds,method,measure,cfg);
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns,scales)plotMeanSTD_1D_EMG_JLDparams(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,scales);
plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns,scales)plotMeanSTD_0D(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,scales);
%freqs = [2:2:20];
%plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns,scales,freqs)plotMeanSTD_2D_JLDparams(mean,std,Measure, Colors, PlotGroups, PlotConds, Ns,scales,freqs);
% -plotMeasure: a handle to a function to plot each measure
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
                  
                    
                    for iM = 1:Nm;%...for each measure...
                        
                        %...get this subject's data into temp...
                        temp=eval(tempCommand{iM});
                        
                        %                         disp(Measures{iM})
                        %                         disp(size(temp))
                        
                        %...finally add this subject's data into a row of PLS
                        %matrix:
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
                        end
                        try
                        PLSmat.(Measures{iM}){iG}(RowIND,:) = temp(:);
                        catch
                            keyboard
                        end
                        
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
            PLSmatStats.(Measures{iM}).mean{iG,iC} = mean( PLSmat.(Measures{iM}){iG}(RowIND,:) );
            PLSmatStats.(Measures{iM}).std{iG,iC} = std( PLSmat.(Measures{iM}){iG}(RowIND,:) );
            
            %...reshape to original sizes of the measures...
            PLSmatStats.(Measures{iM}).mean{iG,iC} = reshape( PLSmatStats.(Measures{iM}).mean{iG,iC}, MeasOutsStruct{iM} );
            PLSmatStats.(Measures{iM}).std{iG,iC} = reshape( PLSmatStats.(Measures{iM}).std{iG,iC}, MeasOutsStruct{iM} );
            
        end
    end%for Nc
end %for Ng
save([PLSmatFolder,filesep,'PLSmat.mat'],'PLSmat','PLSmatStats','Groups','Conds','Subjs','SubjsFiles','Measures','MeasOutsStruct',...
    'Ng','Nc','Ns','Nm','Ntr')
if exist('cfg','var')
    cfgMeasure = cfg;
    save([PLSmatFolder,filesep,'PLSmat.mat'],'cfgMeasure','-append')
end
if exist('Colors','var')
    save([PLSmatFolder,filesep,'PLSmat.mat'],'Colors','-append')
end

%---------------------For Plotting-----------------------------------------
for iM=1:Nm;
    %h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM},Colors, PlotGroups, PlotConds, Ns,cfg, method);
    %h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Ns,Ng,Nc,PlotGroups,PlotConds,method,Measures{iM},cfg);
    h=plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM}, Colors, PlotGroups, PlotConds, Ns,scales);
    %h=plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM}, Colors, PlotGroups, PlotConds, Ns,cfgJLD.scales*0.004,freqs);
    fileName = [PLSmatFolder,filesep,Measures{iM},'_stats_meanSTDerr.fig'];
    saveas(h,fileName);
    close(h);
    clear h;
end
%---------------------For Plotting-----------------------------------------

%---------------------For Plotting-----------------------------------------
for iM=1:Nm;
    %h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM},Colors, PlotGroups, PlotConds, Ns,cfg, method);
    %h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Ns,Ng,Nc,PlotGroups,PlotConds,method,Measures{iM},cfg);
    h=plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM}, Colors, PlotGroups, PlotConds, [],scales);
    %h=plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM}, Colors, PlotGroups, PlotConds, [],cfgJLD.scales*0.004,freqs);
    fileName = [PLSmatFolder,filesep,Measures{iM},'_stats_mean.fig'];
    saveas(h,fileName);
    close(h);
    clear h;
end
%---------------------For Plotting-----------------------------------------

toc

clear all;





%-------------------------B. Run PLS--------------------------------------
%Inputs:
OutputFolder = 'C:\Users\vernooij.c\Documents\MATLAB\FCD\Results\PLS';
%-PLSmatFolder: full path to the folder where PLS matrix is saved
PLSmatFolder = 'C:\Users\vernooij.c\Documents\MATLAB\FCD\Results\PLS';
%-OutputFolder: full path to the folder where PLS output will be stored
%PLSfolder='C:\Users\CoordAgeEEG\Dropbox\DPtoolbox\Statistics\external\plscmd';
PLSfolder='C:\Users\vernooij.c\Documents\Dropbox\DPtoolbox\Statistics\external\plscmd';
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
% load('C:\data\DATA_Aging_Complexity\Results_PSD_JLD\PLS\Aging\Contrast\DesignData.mat')
% load('C:\data\DATA_Aging_Complexity\Results_FC_FCD_JLD\JLD\PLS\RestEO\Aging\Contrast\DesignData.mat')
%  opt.stacked_designdata = DesignData;
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

%Run PLS
addpath(PLSfolder); %add PLS cmd path
disp('Running PLS for measure:...')
for iM = 1:Nm;%...for each measure...
    
    disp(['...',Measures{iM},'...'])
    tic
    %Run PLS
    %result = pls_analysis(datamat_lst, num_subj_lst, num_cond, option)
%    opt.stacked_designdata = DesignData.(Measures{iM});
    PLSres.(Measures{iM}) = pls_analysis(PLSmat.(Measures{iM}), Ns, Nc, opt);
    toc
    
    %Store PLSmat and Msiz(iM)
    PLScfg.(Measures{iM}).PLSmat = PLSmat.(Measures{iM});
    PLSmat = rmfield(PLSmat,Measures{iM});
    PLScfg.(Measures{iM}).MeasOutsStruct = MeasOutsStruct{iM};
    
    
    %Create meaningful ouputs out of this result:
    
%     %A. Percentage of variance explained by each Latent Variable
%     s2 = (PLSres.(Measures{iM}).s).^2; %square of eigenvalues
%     PLSres.(Measures{iM}).var100 = 100*s2/sum(s2);     
            
    
   %B. Brain scores
    Nlv = size(PLSres.(Measures{iM}).s,1);% the number of LVs   
    Ndof = size(PLSres.(Measures{iM}).v,1);%number of degrees of freedom plus 1 (groups x conditions x behavioral measures)
    Nr = size(PLSres.(Measures{iM}).usc,1); % the number of brain scores
    %reshape and normalize brain scores
    uscNorm = PLSres.(Measures{iM}).usc./repmat((PLSres.(Measures{iM}).s).',[Nr,1]);
    %separate brain scores for each group and condition
    for iG = 1:Ng;
        for iC = 1:Nc;
            PLSres.(Measures{iM}).BrainScores{iG,iC} = uscNorm(Nc*sum(Ns(1:iG-1))+(iC-1)*Ns(iG)+[1:Ns(iG)],:);
        end
    end
    
    %B. Brain LVs 
    for iLV=1:Nlv;
            PLSres.(Measures{iM}).BrainLVs{iLV} = reshape(PLSres.(Measures{iM}).u(:,iLV),PLScfg.(Measures{iM}).MeasOutsStruct);
            PLSres.(Measures{iM}).BrainLVsBootsRatios{iLV} = [];
    end 
    %...if there is a bootstrap result (get also brain and task LVs' standard errors)...
    if isfield(PLSres.(Measures{iM}),'boot_result')
 
        for iLV=1:Nlv;
            PLSres.(Measures{iM}).BrainLVsBootsRatios{iLV} = reshape(PLSres.(Measures{iM}).boot_result.compare_u(:,iLV),PLScfg.(Measures{iM}).MeasOutsStruct);
        end
        %reshape the boostrap standard error and error corrected brain scores
        %PLSres.(Measures{iM}).boot_result.compare_u  = reshape(PLSres.(Measures{iM}).boot_result.compare_u,[PLSres.(Measures{iM}).MeasOutsStruct, Nlv]);
        PLSres.(Measures{iM}).boot_result.u_se  = reshape(PLSres.(Measures{iM}).boot_result.u_se,[PLScfg.(Measures{iM}).MeasOutsStruct, Nlv]);
        
        %lower and upper bootstrap standard errors of task LVs v
        PLSres.(Measures{iM}).vlerr = abs( PLSres.(Measures{iM}).boot_result.llusc ./ repmat( (PLSres.(Measures{iM}).s).',[Ndof,1] ) - PLSres.(Measures{iM}).v ); %also normalize with s
        PLSres.(Measures{iM}).vuerr = abs( PLSres.(Measures{iM}).boot_result.ulusc ./ repmat( (PLSres.(Measures{iM}).s).',[Ndof,1] ) - PLSres.(Measures{iM}).v ); %also normalize with s
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




%%-------------------------C. Plotting PLS---------------------------------
scales = unique(round(exp(log(1):(log(1000)-log(1))/100:log(1000))))/100.';
% freqs = [2:2:20];
%Inputs:
OutputFolder = 'C:\Users\vernooij.c\Documents\MATLAB\FCD\Results\PLS';
%-OutputFolder: full path to the folder where PLS output will be stored
alpha = 0.05; %statistical alpha value for significance
%-alpha: alpha value for significance of permutation test of LVs (real 0<scalar<1)
%varMin = 0; %minimum % variance explained by a LV to be plotted
%Method for plotting brain LVs:
%-'BootsRatios': plot bootstrap ratios (possible only if there is a
%                bootstrap test done)
%-'LV': plot the Brain LV
BrainLVplotOpts.method = 'BootsRatios';
%threshold for reliability
BrainLVplotOpts.th = 0; %3 for ~99% confidence
%alpha range for points below threshold for maps of 2D Brain LVs
BrainLVplotOpts.Range = [0 1]; 
Groups={'Young','Old'};%'YC','OC',
Conds={'ID1','ID2','ID3','ID4','ID5'};%;{'R','OnC','OC'}; %'RestEO',
%Colors = {rgb('DodgerBlue'),rgb('ForestGreen'),rgb('Red'), rgb('SteelBlue'),rgb('Teal'),rgb('DarkRed')};
Colors = {rgb('blue'),rgb('green'),rgb('yellow'),rgb('red'),rgb('black')};
%plotBrainLVs = @(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,cfgMeasure)plotUnivComplmap(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,cfgMeasure); 
plotBrainLVs = @(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,scales,Color,LineStyle)plotBrainLVs_1D_EMG_JLDparams(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,scales,Color,LineStyle);
%plotBrainLVs = @(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,scales,freqs)plotBrainLVs_2D_JLDparams(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,scales,freqs);

%varargin={}; %any other inputs you would like to give to the plotting function
%plotBrainLVs = @(BrainLV,Measure,cfgMeasure,varargin)plotBrainLVs_2D(BrainLV,Measure,cfgMeasure,varargin); 
% plotBrainLVs = @(BrainLV,Measure,cfgMeasure,varargin)plotBrainLVs_1D(BrainLV,Measure,cfgMeasure,varargin);
%if you commend this, commend also lines 543 (hb=figure...)-550 (close(hb);)
%-plotBrainLVs: a handle to a function to plot the brain LVs of each measure
%----------------------------End of Inputs---------------------------------

cd(OutputFolder) 
load('PLSres.mat')
%cfgMeasure=PLScfg.cfgMeasure;

if isfield(PLScfg,'Colors')
    Colors = PLScfg.Colors;
end
iX=0;
if PLScfg.opt.method<3
    for iG=1:PLScfg.Ng;
        for iC=1:PLScfg.Nc;
            iX=iX+1;
            xlbl{iX} = [Groups{iG},'-',Conds{iC}];
        end
    end
else
    for iG=1:PLScfg.Ng;
        for iC=1:PLScfg.Nc;
            for iB=1:PLScfg.Nb;
                iX=iX+1;
                xlbl{iX} = [Groups{iG},'-',Conds{iC},'-',Behav{iB}];
            end
        end
    end
end


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
        indLVsign = (thisPLS.perm_result.sprob <= alpha); % & (thisPLS.var100 >= varMin);
        LVsign = LVsign(indLVsign);
        NlvSign = length(LVsign);
    end
    
    
    h=figure('Name',[PLScfg.Measures{iM}, ', task saliences and brain scores']);

    for iLV = 1:NlvSign;
        
        LV  = LVsign(iLV);
             
        figure(h);
        
        subplot(NlvSign,4,4*(iLV-1)+[1 2])
%         if any(PLScfg.opt.method==[1,3]) %if the method is mean-centered...
%             title([PLScfg.Measures{iM},', LV',num2str(LV),', (s = ',num2str(thisPLS.var100(LV)),'%, p = ',num2str(thisPLS.perm_result.sprob(LV)),'): task LVs'])
%         else %if it is non-rotated...
            title([PLScfg.Measures{iM},', LV',num2str(LV),', (s = ',num2str(thisPLS.s(LV)),', p = ',num2str(thisPLS.perm_result.sprob(LV)),'): task LVs'])
%        end
        for iX=1:Ndof;
            hold on;
            bar(iX,thisPLS.v(iX,LV),0.5,'FaceColor',Colors{iX});
            if ( isfield(thisPLS,'vuerr') && isfield(thisPLS,'vlerr') )
                hold on;
                errorbar(iX,thisPLS.v(iX,LV),thisPLS.vlerr(iX,LV),thisPLS.vuerr(iX,LV),'Color','k')
            end
        end
        set(gca,'xlim',[0 Ndof+1],'xtick',[1:Ndof],'xticklabel',xlbl);
        hold off
        
        subplot(NlvSign,4,4*(iLV-1)+[3 4])
        title([PLScfg.Measures{iM},', LV',num2str(LV),': brain scores'])
        hold on;
        iX=0;
        minVal = +inf;
        maxVal =-inf;
        for iG=1:PLScfg.Ng;
            for iC=1:PLScfg.Nc;
                iX=iX+1;
                temp = double(thisPLS.BrainScores{iG,iC}(:,LV));
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
        hold off
         
        hb=figure('Name',[PLScfg.Measures{iM},', LV',num2str(LV),': brain latent variable']);
        %         if any(PLScfg.opt.method==[1,3]) %if the method is mean-centered...
        %title([PLScfg.Measures{iM},', LV',num2str(LV),', (s = ',num2str(thisPLS.var100(LV)),'%, p = ',num2str(thisPLS.perm_result.sprob(LV)),'): brain LVs'])
        %else
        title([PLScfg.Measures{iM},', LV',num2str(LV),', (s = ',num2str(thisPLS.s(LV)),', p = ',num2str(thisPLS.perm_result.sprob(LV)),'): brain LVs'])
        %end        hold on;
        %plotBrainLVs(thisPLS.BrainLVs{LV},thisPLS.BrainLVsBootsRatios{LV},BrainLVplotOpts,PLScfg.Measures{iM},cfgMeasure);
        plotBrainLVs(thisPLS.BrainLVs{LV},thisPLS.BrainLVsBootsRatios{LV},BrainLVplotOpts,PLScfg.Measures{iM},scales,'k','-');
        %plotBrainLVs(thisPLS.BrainLVs{LV},thisPLS.BrainLVsBootsRatios{LV},BrainLVplotOpts,PLScfg.Measures{iM},scales,freqs);
        hold off;
        fileNameB = [PLScfg.Measures{iM},'_LV',num2str(LV),'_LV'];
        saveas(hb,[fileNameB,'.fig']);
        close(hb);
        
    end
    fileName = [PLScfg.Measures{iM},'_TaskLVs_Scores'];
    saveas(h,[fileName,'.fig']);
%     close(h);
end

