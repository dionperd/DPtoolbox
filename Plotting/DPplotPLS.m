%% %-------------------------C. Plotting PLS---------------------------------
%Inputs:
OutputFolder=pwd;%'C:\Users\dionperd\Desktop\wetransfer-613a7c';
FileName = 'PLSres.mat';%'PLSres_al_to_re_500.mat'; %Default 'PLSres.mat'
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
BrainLVplotOpts.th = 0; %3 stands for ~99% confidence
%alpha range for points below threshold for maps of 2D Brain LVs
BrainLVplotOpts.Range = [0 0.5]; 
%-alpha: alpha value for significance of permutation test of LVs (real 0<scalar<1)

Groups = {};%{'G1','G2'};
Conds = {};%{'C1','C2','C3'};%
Behav = {};
cfgMeasure = [];
LineColor = 'b';
LineStyle = '-';
plotBrainLVs = @(BrainLV,BrainLVbootsRatio,BrainLVplotOpts,Measure,LineColor,LineStyle,cfgMeasure)plotBrainLVs_2D(BrainLV,BrainLVbootsRatio,BrainLVplotOpts,Measure,LineColor,LineStyle,cfgMeasure); 
%plotBrainLVs = @(BrainLV,BrainLVbootsRatio,BrainLVplotOpts,Measure,LineColor,LineStyle,cfgMeasure)plotBrainLVs_1D(BrainLV,BrainLVbootsRatio,BrainLVplotOpts,Measure,LineColor,LineStyle,cfgMeasure); 
%-plotBrainLVs: a handle to a function to plot the brain LVs of each measure
%----------------------------End of Inputs---------------------------------

cd(OutputFolder) 
load(FileName)
Ng = PLScfg.Ng;
Nc = PLScfg.Nc;
Nb = PLScfg.Nb;

if isempty(Groups)
    if isfield(PLScfg,'Groups')
        Groups  = PLScfg.Groups;
    else
        for iG=1:Ng;
            Groups{iG} =['G',num2str(iG)];
        end
    end
end
if isempty(Conds)
    if isfield(PLScfg,'Conds')
        Conds  = PLScfg.Conds;
    else
        for iC=1:Nc;
            Conds{iC} =['C',num2str(iC)];
        end
    end
end
if isempty(Behav)
    if isfield(PLScfg,'Behav')
        Behav  = PLScfg.Behav;
    else
        for iB=1:Nb;
            Behav{iB} =['B',num2str(iB)];
        end
    end
end

if isfield(PLScfg,'Colors')
    Colors = PLScfg.Colors;
else
    Colors = {rgb('DodgerBlue'),rgb('ForestGreen'),rgb('Red'), rgb('DarkBlue'),rgb('Teal'),rgb('DarkRed')};
end
if isfield(PLScfg,'cfgMeasure')
    cfgMeasure = PLScfg.cfgMeasure;   
end

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


for iM = 1:PLScfg.Nm;%...for each measure...
    
    disp('Plotting PLS results for measure:...')
    disp(PLScfg.Measures{iM})
    
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
    
    
    h=figure('Visible','off','Name',[PLScfg.Measures{iM}, ', task saliences and brain scores']);

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
            if isfield(thisPLS,'p')
                title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, ',num2str(thisPLS.var100(LV)),'%, p < ',num2str(thisPLS.p(LV)),'):'];['task LVs and standardized brain scores'' distribution']})
            else
                title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, ',num2str(thisPLS.var100(LV)),'%):'];['task LVs and standardized brain scores'' distribution']})
            end
        else %if it is non-rotated...
            if isfield(thisPLS,'p')
                title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, s = ',num2str(thisPLS.s(LV)),', p < ',num2str(thisPLS.p(LV)),'):'];['task LVs and standardized brain scores'' distribution']})
            else
                title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, s = ',num2str(thisPLS.s(LV)),'):'];['task LVs and standardized brain scores'' distribution']})                
            end
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
            if isfield(thisPLS,'p')
                title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, ',num2str(thisPLS.var100(LV)),'%, p < ',num2str(thisPLS.p(LV)),'):'];['brain scores']})
            else
                title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, ',num2str(thisPLS.var100(LV)),'%):'];['brain scores']})                
            end
        else %if it is non-rotated...
            if isfield(thisPLS,'p')
                title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, s = ',num2str(thisPLS.s(LV)),', p < ',num2str(thisPLS.p(LV)),'):'];['brain scores']})
            else
                title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, s = ',num2str(thisPLS.s(LV)),'):'];['brain scores']})                
            end
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
            if isfield(thisPLS,'p')
                title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, ',num2str(thisPLS.var100(LV)),'%, p < ',num2str(thisPLS.p(LV)),'):'];['brain scores per subject']})
            else
                title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, ',num2str(thisPLS.var100(LV)),'%):'];['brain scores per subject']})                
            end
        else %if it is non-rotated...
            if isfield(thisPLS,'p')
                title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, s = ',num2str(thisPLS.s(LV)),', p < ',num2str(thisPLS.p(LV)),'):'];['brain scores per subject']})
            else
                title({[PLScfg.Measures{iM},' (LV_{',num2str(LV),'}, s = ',num2str(thisPLS.s(LV)),'):'];['brain scores per subject']})                
            end
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
         
        hb=figure('Visible','on','Name',[PLScfg.Measures{iM},', LV',num2str(LV),': brain latent variable'],'units','normalized','outerposition',[0 0 1 1]);
        if any(PLScfg.opt.method==[1,3]) %if the method is mean-centered...
            if isfield(thisPLS,'p')
                title([PLScfg.Measures{iM},', brain LV_{',num2str(LV),'} (',num2str(thisPLS.var100(LV)),'%, p < ',num2str(thisPLS.p(LV)),')'])
            else
                title([PLScfg.Measures{iM},', brain LV_{',num2str(LV),'} (',num2str(thisPLS.var100(LV)),'%)'])
            end
        else
            if isfield(thisPLS,'p')
                title([PLScfg.Measures{iM},', brain LV_{',num2str(LV),'} (s = ',num2str(thisPLS.s(LV)),', p < ',num2str(thisPLS.p(LV)),')'])
            else
                title([PLScfg.Measures{iM},', brain LV_{',num2str(LV),'} (s = ',num2str(thisPLS.s(LV)),')'])
            end
        end
        hold on;
        plotBrainLVs(sgn*thisPLS.BrainLVs{LV},sgn*thisPLS.BrainLVsBootsRatios{LV},BrainLVplotOpts,PLScfg.Measures{iM},LineColor,LineStyle,cfgMeasure);
        hold off;
        fileNameB = [PLScfg.Measures{iM},'_LV',num2str(LV),'_BrainLV'];
        saveas(hb,[fileNameB,'.fig']);
        %close(hb);
        
    end
    fileName = [PLScfg.Measures{iM},'_TaskLVs_BrainScores'];
    saveas(h,[fileName,'.fig']);
    %close(h);
end