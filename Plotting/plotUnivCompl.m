function h=plotUnivCompl(LV,BS,measure,cfg,Color,LineStyle,opts)

if nargin<5
    LineStyle = '-';
end
if nargin<4
    Color='b';
end

if strcmpi(measure,'MSE')
    x=cfg.MSE.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'MSE';
    h=plotElectrodes(x,LV,BS,xlbl,ylbl,Color,LineStyle,opts);
    
elseif strcmpi(measure,'CMSE')
    x=cfg.CMSE.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'CMSE';
    h=plotElectrodes(x,LV,BS,xlbl,ylbl,Color,LineStyle,opts);
    
elseif strcmpi(measure,'STD')
    x=cfg.MSE.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'STD';
    h=plotElectrodes(x,LV,BS,xlbl,ylbl,Color,LineStyle,opts);
    
elseif strcmpi(measure,'logF')
    x=log(cfg.DFA.scales/cfg.fs);
    xlbl = 'log(scales) (log(sec))';
    ylbl = 'logF';
    h=plotElectrodes(x,LV,BS,xlbl,ylbl,Color,LineStyle,opts);
    
elseif strcmpi(measure,'H')
    x=cfg.DFA.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'H';
    h=plotElectrodes(x,LV,BS,xlbl,ylbl,Color,LineStyle,opts);    
    
elseif strcmpi(measure,'logV')
    x=log(cfg.VGR.scales/cfg.fs);
    xlbl = 'log(scales) (log(sec))';
    ylbl = 'log(VAR)';
    h=plotElectrodes(x,LV,BS,xlbl,ylbl,Color,LineStyle,opts);
    
elseif strcmpi(measure,'logP')
    x=cfg.SP.logf;
    xlbl = 'log(f) (log(Hz))';
    ylbl = 'log(POW)';
    plotElectrodes(x,LV,BS,xlbl,ylbl,Color,LineStyle,opts);
    
elseif strcmpi(measure,'P')
    x= cfg.SP.f(2:end);
    xlbl = 'f (Hz)';
    ylbl = 'POW';
    h=plotElectrodes(x,LV,BS,xlbl,ylbl,Color,LineStyle,opts);

elseif strcmpi(measure,'PS')||strcmpi(measure,'VS')||strcmpi(measure,'DOF')
    h=topoplot(LV,BS,opts);
    
end


function h=plotElectrodes(x,LV,BS,xlbl,ylbl,Color,LineStyle,opts)

%19 channels, Cz=10:
Ch=[1,3,6,8,10,12,14,22,24,26,28,30,40,42,44,46,48,56,58]; %electrodes
ChLab = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8','P7','P3','Pz','P4','P8','O1','O2'};

%58 channels, Cz=26:
% Ch=[1:58]; %electrodes
% ChLab = { 'Fp1', 'Fpz', 'Fp2', 'AF3', 'AF4', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4',...
%           'F6', 'F8', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'T7', 'C5', 'C3',...
%           'C1', 'Cz', 'C2', 'C4', 'C6', 'T8', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2',...
%           'CP4', 'CP6', 'TP8', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8',...
%           'PO7', 'PO3', 'PO1', 'POz', 'PO2', 'PO4', 'PO8', 'O1', 'Oz', 'O2'};
  
if isempty(BS)
    y = LV;
else
    inds = abs(BS)>opts.th;
    if strcmpi(opts.method,'BootsRatios')
        y=BS;
    else
        y=LV;
    end
end
y = y(:,Ch);
inds = inds(:,Ch);
Nc=length(Ch);
SubPlots=[2,4,6:20,22,24];

minVal = min(y(:));
maxVal = max(y(:));
for iC=1:Nc;
    subplot(5,5,SubPlots(iC));
    hold on;
    plot(x,y(:,Ch(iC)),'color',Color,'LineStyle',LineStyle,'linewidth',2);
    if strcmpi(opts.method,'BootsRatios')
        plot(x,y(:,Ch(iC)),'color',Color,'LineStyle',LineStyle,'linewidth',2);
        plot(x,opts.th*ones(1,length(x)),'--','color',Color,'LineStyle',LineStyle,'linewidth',0.5);
        plot(x,-opts.th*ones(1,length(x)),'--','color',Color,'LineStyle',LineStyle,'linewidth',0.5);
    else
        plot(x( inds(:,Ch(iC)) ), y( inds(:,Ch(iC)), Ch(iC) ),'*','markersize',2,'markerfacecolor',Color,'markeredgecolor',Color);
    end
    set(gca,'ylim',[minVal,maxVal]);
    grid on;
    title(ChLab{iC})
    xlabel(xlbl);ylabel(ylbl);
   
end
h=gcf;


function  h=topoplot(LV,BS,opts)

% %Ch=[1,3,6,8,10,12,14,22,24,26,28,30,40,42,44,46,48,56,58]; %electrodes
% 
% chx = [ 53.0000   53.0000  107.0000  121.0000  125.0000  122.0000  114.0000  204.0000  199.0000...
%         198.0000  200.0000  206.0000  288.0000  276.0000  271.0000  276.0000  288.0000  341.0000  342.0000].';
%     
% chy = [ 148.0000  239.0000   71.0000  133.0000  192.0000  255.0000  313.0000   51.0000  119.0000...
%         191.0000  259.0000  332.0000   75.0000  128.0000  190.0000  253.0000  305.0000  147.0000  231.0000].';
% LV = LV(Ch);
% Ch=[1:19]; %electrodes

chx = [ 50.9956,   57.9605,   54.4781,   86.9810,   86.9810,  107.8757,  114.8406,  120.6447,  124.1272,  126.4488, 125.2880,  121.8056,  118.3231,  111.3582,  155.4693,  158.9518,  163.5950,  162.4342,  162.4342,  161.2734, 156.6301,  200.7412,  200.7412,  198.4196,  197.2588,  197.2588,  196.0980,  198.4196,  198.4196,  207.7061, 246.0132,  246.0132,  239.0482,  236.7266,  235.5658,  237.8874,  241.3699,  246.0132,  248.3348,  283.1594, 279.6769,  273.8728,  271.5512,  269.2295,  272.7120,  273.8728,  284.3202, 287.8026,  316.8231,  312.1798, 301.7325,  306.3757,  304.0541, 313.3406,  321.4664,  340.0395,  341.2003,  342.3611].';

chy = [  146.3710  193.3387  237.5968,  135.5323  250.2419   70.5000  100.3065  131.0161  163.5323  190.6290 223.1452  254.7581  283.6613  316.1774   88.5645  123.7903  157.2097  190.6290  224.9516  258.3710 294.5000   50.6290   85.8548  119.2742  153.5968  191.5323  227.6613  261.0806  297.2097  334.2419 52.4355   92.1774  120.1774  152.6935  190.6290  223.1452  255.6613  291.7903  323.4032   75.0161 104.8226  125.5968  158.1129  190.6290  220.4355  252.0484  278.2419  305.3387  103.0161  136.4355 166.2419  187.0161  214.1129  243.9194  272.8226  148.1774  187.9194  229.4677 ]';

%Ch=[1:58]; %electrodes

if isempty(BS) 
    x=LV(:);
    minVal = min(x);
    maxVal = max(x);
elseif strcmpi(opts.method,'BootsRatios')
    x=BS(:);
    minVal = -5;
    maxVal = 5;
else 
    x=LV(:);
    minVal = min(x);
    maxVal = max(x);
end
% oldFolder=cd('C:\Users\dionperd\Dropbox\Dionysis\DPtoolbox\Plotting\eegplot');
%oldFolder=cd('C:\Users\CoordAgeEEG\Dropbox\DPtoolbox\Plotting\eegplot');
h=eegplot(double(x),[chx, chy],1,1,'nearest',[],minVal,maxVal);

%cd(oldFolder)