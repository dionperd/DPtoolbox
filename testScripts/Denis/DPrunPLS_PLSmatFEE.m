%function [PLSres, PLScfg] = DPrunPLStask( PLSfolder, RESfolder,  RESstruct, Conds, Groups, Measures, MeasOutsStruct, Ntr, Ns, Subjs)
%DPrunPLS This function calculates PLS analysis for several measures on a
%potentially multigroup and multisubject dataset

%   Outputs:
% -PLSres: structure with PLS results for each measure
% -PLScfg: structure with PLS configuration (PLSmat included)  for each measure

%---------------------A. Stack PLS matrix----------------------------------
%   Inputs:
RESfolder='/Users/dionysos/Dropbox/Dionysis/postdocberlin/Publications/BrainSynchronyFEEs/brain_data/DvsU';
%-RESfolder: full path to the folder of results .mat files
PLSmatFolder = '/Users/dionysos/Dropbox/Dionysis/postdocberlin/Publications/BrainSynchronyFEEs/brain_data/DvsU';
%-PLSmatFolder: full path to the folder where PLS matrix will be saved
PLSmat = plsDvsU;
%clear plsDvsS;
%-RESstruct: name of the result structure (string)
Conds = PLSmat.cfg.Conds;%
%Conds={'Cond1','Cond2','Cond3'};%
%-Conds: cellstring of the names of conditions
Groups = PLSmat.cfg.Groups;%
%Groups={'Group1','Group2'};
%-Groups: cellstring of the names of Groups 
%(a string common to all files, e.g., the extension, if there is only 1 group)
Ns = PLSmat.cfg.Ns;
%Ns= [31 28 40]; 
%-Ns: number of subjects of each group (vector of positive integers)
%load('C:\Users\CoordAgeEEG\Dropbox\Results\Subjects.mat')
%Subjs =  PLSmat.cfg.Subjs;%
%-Subj: cell of cellstrings of the names of subjects for each group in case we want to use a
%       specific subset of them, otherwise Subjs={}
Measures = {'WP','PLI'};
%-Measures: a cellstring of the names of measures to be calculated
MeasOutsStruct = {[38 16 60], [38 16 60]};
%-MeasOutsStruct: a cell of the same size as Measures with vectors
%                 corresponding to the size of the specific measure's 
%                 result PER TRIAL


%---------------------For Plotting-----------------------------------------
%Inputs for plotting:
Colors = {rgb('DodgerBlue'),rgb('ForestGreen'), rgb('DarkBlue'),rgb('Teal')};
varargin={}; %any other inputs you would like to give to the plotting function
plotMeasure = @(mean,std,Measure, Colors, Groups, Conds, Ns,varargin)plotMeanSTD_3D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,varargin)plotMeanSTD_1D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,varargin)plotMeanSTD_0D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
% % -plotMeasure: a handle to a function to plot each measure
%--------------------For Plotting------------------------------------------
%----------------------------End of Inputs---------------------------------




Nc = max(numel(Conds),1);
Ng = max(numel(Groups),1);
Nm =  numel(Measures);


Nr = Ns*Nc; %number of rows of PLSmat
%Initialize PLS matrices for each (sub)measure
for iM = 1:Nm;%For each measure...
    PLSmatStats.(Measures{iM}).mean = cell(Ng,Nc);
    PLSmatStats.(Measures{iM}).std = cell(Ng,Nc);
end



tic
for iG = 1:Ng; %For each group...
    disp('Group=')
    disp(Groups{iG})
    for iC = 1:Nc; %...for each condition...
        disp('Condition=')
        disp(Conds{iC})
        
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

save([PLSmatFolder,filesep,'PLSmat.mat'],'PLSmat','PLSmatStats','Groups','Conds','Measures','MeasOutsStruct',...
    'Ng','Nc','Ns','Nm')
if exist('cfg','var')
    cfgMeasure = cfg;
    save('PLSmat.mat','cfgMeasure','-append')
end
if exist('Colors','var')
    save('PLSmat.mat','Colors','-append')
end

% %---------------------For Plotting-----------------------------------------
% for iM=1:Nm;
%     h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM}, Colors, Groups, Conds, []);
%     fileName = [PLSmatFolder, filesep, Measures{iM},'_mean.fig'];
%     saveas(h,fileName);
%     close(h);
%     clear h;
% end
% %---------------------For Plotting-----------------------------------------
% 
% %---------------------For Plotting-----------------------------------------
% for iM=1:Nm;
%     h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM}, Colors, Groups, Conds, Ns);
%     fileName = [PLSmatFolder, filesep, Measures{iM},'_meanSTDerr.fig'];
%     saveas(h,fileName);
%     close(h);
%     clear h;
% end
% %---------------------For Plotting-----------------------------------------


toc

clear all;





%%-------------------------B. Run PLS--------------------------------------
%Inputs:
PLSmatFolder ='/Users/dionysos/Dropbox/Dionysis/postdocberlin/Publications/BrainSynchronyFEEs/brain_data/DvsU';
%-PLSmatFolder: full path to the folder where PLS matrix is saved
OutputFolder='/Users/dionysos/Dropbox/Dionysis/postdocberlin/Publications/BrainSynchronyFEEs/brain_data/DvsU';
%-OutputFolder: full path to the folder where PLS output will be stored
PLSfolder='C:\Users\CoordAgeEEG\Dropbox\DPtoolbox\Statistics\external\plscmd';
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

%Run PLS
addpath(PLSfolder); %add PLS cmd path
disp('Running PLS for measure:...')
for iM = 1:Nm;%...for each measure...
    
    disp(['...',Measures{iM},'...'])
    tic
    %Run PLS
    %result = pls_analysis(datamat_lst, num_subj_lst, num_cond, option)
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
%PLScfg.Ntr = Ntr;
PLScfg.Nlv =Nlv;
PLScfg.Ndof = Ndof;
PLScfg.Conds = Conds;
PLScfg.Groups = Groups;
if opt.method==3
    PLScfg.Behav = Behav;
end
% PLScfg.Subjs = Subjs;
% PLScfg.SubjsFiles = SubjsFiles;
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





% %%-------------------------C. Plotting PLS---------------------------------
% %Inputs:
% OutputFolder='/Users/dionysos/Dropbox/Dionysis/postdocberlin/Publications/BrainSynchronyFEEs/brain_data/DvsS';
% %-OutputFolder: full path to the folder where PLS output will be stored
% alpha = 0.1; %statistical alpha value for significance
% %-alpha: alpha value for significance of permutation test of LVs (real 0<scalar<1)
% %varMin = 1; %minimum % variance explained by a LV to be plotted
% %Method for plotting brain LVs:
% %-'BootsRatios': plot bootstrap ratios (possible only if there is a
% %                bootstrap test done)
% %-'LV': plot the Brain LV
% BrainLVplotOpts.method = 'BootsRatios';
% %threshold for reliability
% BrainLVplotOpts.th = 0; %~99% confidence
% %alpha range for points below threshold for maps of 2D Brain LVs
% BrainLVplotOpts.Range = [0 1]; 
% %-alpha: alpha value for significance of permutation test of LVs (real 0<scalar<1)
% Groups = {''};
% %Groups = {'OwithAlpha','OnoAlpha'};
% Conds={'aD','hD','aS','hS'};%
% %Behav = {'Age'};
% 
% Colors = {rgb('DodgerBlue'),rgb('ForestGreen'), rgb('DarkBlue'),rgb('Teal')};
% cfgMeasure = PLScfg.cfgMeasure;
% plotBrainLVs = @(BrainLV,BrainLVbootsRatio,BrainLVplotOpts,Measure,cfgMeasure,varargin)plotBrainLVs_2D(BrainLV,BrainLVbootsRatio,BrainLVplotOpts,Measure,cfgMeasure); 
% % plotBrainLVs = @(BrainLV,Measure,cfgMeasure)plotBrainLVs_1D(BrainLV,BrainLVbootsRatio,BrainLVplotOpts,Measure,cfgMeasure,varargin); 
% %-plotBrainLVs: a handle to a function to plot the brain LVs of each measure
% %----------------------------End of Inputs---------------------------------
% 
% cd(OutputFolder) 
% load('PLSres.mat')
% 
% if isfield(PLScfg,'Colors')
%     Colors = PLScfg.Colors;
% end
% iX=0;
% if PLScfg.opt.method<3
%     for iG=1:PLScfg.Ng;
%         for iC=1:PLScfg.Nc;
%             iX=iX+1;
%             xlbl{iX} = [Groups{iG},'-',Conds{iC}];
%         end
%     end
% else
%     for iG=1:PLScfg.Ng;
%         for iC=1:PLScfg.Nc;
%             for iB=1:PLScfg.Nb;
%                 iX=iX+1;
%                 xlbl{iX} = [Groups{iG},'-',Conds{iC},'-',Behav{iB}];
%             end
%         end
%     end
% end
% 
% Nlv = PLScfg.Nlv;
% Ndof = PLScfg.Ndof;
% disp('Plotting PLS results for measure:...')
% for iM = 1:PLScfg.Nm;%...for each measure...
%     
%     thisPLS = PLSres.(PLScfg.Measures{iM});
% 
%     NlvSign = Nlv;
%     LVsign = 1:NlvSign;
%     if isfield(thisPLS,'perm_result') 
%         %Choose the LVs that have p<= alpha and explain more than varMin %
%         %of variance:
%         indLVsign = (thisPLS.perm_result.sprob <= alpha);% & (thisPLS.var100 >= varMin);
%         LVsign = LVsign(indLVsign);
%         NlvSign = length(LVsign);
%     end
%     
%     
%     h=figure('Name',[PLScfg.Measures{iM}, ', task saliences and brain scores']);
% 
%     for iLV = 1:NlvSign;
%         
%         LV  = LVsign(iLV);
%         
%         
%         figure(h);
%         
%         subplot(NlvSign,4,4*(iLV-1)+[1 2])
% %         if any(PLScfg.opt.method==[1,3]) %if the method is mean-centered...
% %             title([PLScfg.Measures{iM},', LV',num2str(LV),', (s = ',num2str(thisPLS.var100(LV)),'%, p = ',num2str(thisPLS.perm_result.sprob(LV)),'): task LVs'])
% %         else %if it is non-rotated...
%             title([PLScfg.Measures{iM},', LV',num2str(LV),', (s = ',num2str(thisPLS.s(LV)),', p = ',num2str(thisPLS.perm_result.sprob(LV)),'): task LVs'])
% %        end
%         for iX=1:Ndof;
%             hold on;
%             bar(iX,thisPLS.v(iX,LV),0.5,'FaceColor',Colors{iX});
%             if ( isfield(thisPLS,'vuerr') && isfield(thisPLS,'vlerr') )
%                 hold on;
%                 errorbar(iX,thisPLS.v(iX,LV),thisPLS.vlerr(iX,LV),thisPLS.vuerr(iX,LV),'Color','k')
%             end
%         end
%         set(gca,'xlim',[0 Ndof+1],'xtick',[1:Ndof],'xticklabel',xlbl);
%         hold off
%         
%         subplot(NlvSign,4,4*(iLV-1)+[3 4])
%         title([PLScfg.Measures{iM},', LV',num2str(LV),': brain scores'])
%         hold on;
%         iX=0;
%         minVal = +inf;
%         maxVal =-inf;
%         for iG=1:PLScfg.Ng;
%             for iC=1:PLScfg.Nc;
%                 iX=iX+1;
%                 temp = double(thisPLS.BrainScores{iG,iC}(:,LV));
%                 minTemp=min(temp);
%                 if minTemp<minVal
%                     minVal=minTemp;
%                 end
%                 maxTemp=max(temp);
%                 if maxTemp>maxVal
%                     maxVal=maxTemp;
%                 end
%                 for iS = 1:PLScfg.Ns(iG);
%                     text(iX,temp(iS),num2str(iS),'Color',Colors{iX});
%                 end
%             end
%         end
%         set(gca,'xlim',[0 Ndof+1],'ylim',[minVal maxVal],'xtick',[1:Ndof],'xticklabel',xlbl);
%         hold off
%          
%         hb=figure('Name',[PLScfg.Measures{iM},', LV',num2str(LV),': brain latent variable']);
%         %         if any(PLScfg.opt.method==[1,3]) %if the method is mean-centered...
%         %title([PLScfg.Measures{iM},', LV',num2str(LV),', (s = ',num2str(thisPLS.var100(LV)),'%, p = ',num2str(thisPLS.perm_result.sprob(LV)),'): brain LVs'])
%         %else
%         title([PLScfg.Measures{iM},', LV',num2str(LV),', (s = ',num2str(thisPLS.s(LV)),', p = ',num2str(thisPLS.perm_result.sprob(LV)),'): brain LVs'])
%         %end
%         hold on;
%         plotBrainLVs(thisPLS.BrainLVs{LV},thisPLS.BrainLVsBootsRatios{LV},BrainLVplotOpts,PLScfg.Measures{iM},cfgMeasure);
%         hold off;
%         fileNameB = [PLScfg.Measures{iM},'_LV',num2str(LV),'_BrainLV'];
%         saveas(hb,[fileNameB,'.fig']);
%         close(hb);
%         
%     end
%     fileName = [PLScfg.Measures{iM},'_TaskLVs_BrainScores'];
%     saveas(h,[fileName,'.fig']);
%     close(h);
% end
