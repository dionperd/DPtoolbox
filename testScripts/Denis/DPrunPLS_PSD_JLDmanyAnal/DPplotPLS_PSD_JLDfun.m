function DPplotPLS_PSD_JLDfun(PLSfolder)


%% %-------------------------C. Plotting PLS---------------------------------
%Inputs:
OutputFolder= PLSfolder;
%-OutputFolder: full path to the folder where PLS output will be stored
alpha = 1; %statistical alpha value for significance
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


Conds={'REC','UOT','AOT'};%'RestEC','Odd_nC','Odd_C'
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
plotBrainLVs = @(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,cfgMeasure)plotBrainLVs_2D_JLDparams(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,cfgMeasure);
%plotBrainLVs = @(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,cfgMeasure)plotBrainLVs_1D_JLDparams_freqs(BrainLV,BrainLVBootsRatios,BrainLVplotOpts,Measure,cfgMeasure);
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
Nm = length(Measures);
Nlv = PLScfg.Nlv;
Ndof = PLScfg.Ndof;
disp('Plotting PLS results for measure:...')

h=figure('Name',['Lifespan EEG synchronization dynamics'],'units','normalized','outerposition',[0 0 1 1]); %'Visible','off',
%,'units','normalized','outerposition',[0 0 1 1]);

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
    
    

    for iLV = 1:1;
        
        LV  = LVsign(iLV);
        
        
        if thisPLS.v(1,iLV) > 0
            sgn = -1;
        else
            sgn = 1;
        end
        
        figure(h);
        
        subplot(3,Nm,iM)
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
        
%         subplot(3,2,3)
%         title('brain scores')
%         for iX=1:Ndof;
%             hold on;
%             if isfield(PLSres.(Measures{iM}),'boot_result')
%                 thisM = (sgn*thisPLS.BrainScoresBootsMeanCntrMU{iLV}(iX)+sgn*thisPLS.BrainScoresGrandMean(iLV));
%                 bar(iX,thisM,0.5,'FaceColor',Colors{iX});
%                 errorbar(iX,thisM,...
%                     min([sgn*thisPLS.BrainScoresBootsLL{iLV}(iX),sgn*thisPLS.BrainScoresBootsUL{iLV}(iX)]),...
%                     max([sgn*thisPLS.BrainScoresBootsLL{iLV}(iX),sgn*thisPLS.BrainScoresBootsUL{iLV}(iX)]),'k')
%             else
%                 bar(iX,sgn*thisPLS.BrainScoresMU{iLV}(iX),0.5,'FaceColor',Colors{iX});
%                 errorbar(iX,sgn*thisPLS.BrainScoresMU{iLV}(iX),...
%                     sgn*thisPLS.BrainScoresMU{iLV}(iX) - thisPLS.BrainScoresSTDERR{iLV}(iX),...
%                     sgn*thisPLS.(Measures{iM}).BrainScoresMU{iLV}(iX) + thisPLS.BrainScoresSTDERR{iLV}(iX),'k')
%             end
%         end
%         set(gca,'xlim',[0 Ndof+1],'xtick',[1:Ndof],'xticklabel',xlbl);
%         for iX = 1:Ndof;
%             if numel(xlbl{iX})>5
%                 xticklabel_rotate([],45)
%                 break;
%             end
%         end
%         hold off
        
        subplot(3,Nm,Nm+iM)
        title('brain scores per subject')
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
        
        
        subplot(3,Nm,2*Nm+iM)
        title(Measures{iM})
        hold on;
        plotBrainLVs(sgn*thisPLS.BrainLVs{LV},sgn*thisPLS.BrainLVsBootsRatios{LV},BrainLVplotOpts,Measures{iM},cfgMeasure);
        hold off;
        
    
    end

end

fileName = 'PLSres';
saveas(h,[fileName,'.fig']);
close(h);