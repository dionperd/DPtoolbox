%-------------------------B. Run PLS--------------------------------------
%Inputs:
PLSmatFolder =pwd;%'full path to your PLS matrix folder';
%-PLSmatFolder: full path to the folder where PLS matrix is saved
OutputFolder=pwd;%'full path to your desired PLS result folder';
if ~exist(OutputFolder,'dir')
    mkdir(OutputFolder)
end
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
    PLSres.(Measures{iM}) = DPpls_analysis(PLSmat.(Measures{iM}), Ns, Nc, opt);
    toc
    
    %Store PLSmat and Msiz(iM)
    PLScfg.(Measures{iM}).PLSmat = PLSmat.(Measures{iM});
    PLSmat = rmfield(PLSmat,Measures{iM});
    PLScfg.(Measures{iM}).MeasOutsStruct = MeasOutsStruct{iM};
    
    
    %Create meaningful ouputs out of this result:
    
% Now this output is returned by DPpls_analysis by default for mean
% centered PLS
%     %A. Percentage of variance explained by each Latent Variable
%     s2 = (PLSres.(Measures{iM}).s).^2; %square of eigenvalues
%     PLSres.(Measures{iM}).var100 = 100*s2/sum(s2);  

   %Calculate p value 
   if isfield(PLSres.(Measures{iM}),'perm_result')
        PLSres.(Measures{iM}).p = max(ceil(PLSres.(Measures{iM}).perm_result.sprob*opt.num_perm)/opt.num_perm,1/opt.num_perm); 
   end
   
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





