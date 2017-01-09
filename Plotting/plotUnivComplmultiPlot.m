function h=plotUnivComplmultiPlot(LV,measure,cfg,Color,LineStyle,sbplt)


if nargin<6
    sbplt=[1 1 1];
end
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
    h=plotElectrodes(x,LV,xlbl,ylbl,Color,LineStyle,sbplt);
    
elseif strcmpi(measure,'CMSE')
    x=cfg.CMSE.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'CMSE';
    h=plotElectrodes(x,LV,xlbl,ylbl,Color,LineStyle,sbplt);
    
elseif strcmpi(measure,'STD')
    x=cfg.MSE.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'STD';
    h=plotElectrodes(x,LV,xlbl,ylbl,Color,LineStyle,sbplt);
    
elseif strcmpi(measure,'logF')
    x=log(cfg.DFA.scales/cfg.fs);
    xlbl = 'log(scales) (log(sec))';
    ylbl = 'logF';
    h=plotElectrodes(x,LV,xlbl,ylbl,Color,LineStyle,sbplt);
    
elseif strcmpi(measure,'H')
    x=cfg.DFA.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'H';
    h=plotElectrodes(x,LV,xlbl,ylbl,Color,LineStyle,sbplt);
    
elseif strcmpi(measure,'logV')
    x=log(cfg.VGR.scales/cfg.fs);
    xlbl = 'log(scales) (log(sec))';
    ylbl = 'log(VAR)';
    h=plotElectrodes(x,LV,xlbl,ylbl,Color,LineStyle,sbplt);
    
elseif strcmpi(measure,'logP')
    x=cfg.SP.logf;
    xlbl = 'log(f) (log(Hz))';
    ylbl = 'log(POW)';
    h=plotElectrodes(x,LV,xlbl,ylbl,Color,LineStyle,sbplt);
    
elseif strcmpi(measure,'P')
    x=cfg.SP.f(2:end);
    xlbl = 'f (Hz)';
    ylbl = 'POW';
    h=plotElectrodes(x,LV,xlbl,ylbl,Color,LineStyle,sbplt);
    
    
elseif strcmpi(measure,'PS')||strcmpi(measure,'VS')||strcmpi(measure,'DOF')||strcmpi(measure,'mMSE')||strcmpi(measure,'mMSEn')
    h=topoplot(LV);
    
end


function h=plotElectrodes(x,LV,xlbl,ylbl,Color,LineStyle,sbplt)

%19 channels, Cz=10:
%Ch=[1,3,6,8,10,12,14,22,24,26,28,30,40,42,44,46,48,56,58]; %electrodes
%ChLab = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8','P7','P3','Pz','P4','P8','O1','O2'};

%58 channels, Cz=26:
% Ch=[1:58]; %electrodes
ChLab = { 'Fp1', 'Fpz', 'Fp2', 'AF3', 'AF4', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4',...
          'F6', 'F8', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'T7', 'C5', 'C3',...
          'C1', 'Cz', 'C2', 'C4', 'C6', 'T8', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2',...
          'CP4', 'CP6', 'TP8', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8',...
          'PO7', 'PO3', 'PO1', 'POz', 'PO2', 'PO4', 'PO8', 'O1', 'Oz', 'O2'};
      
% Nc=length(Ch);

iC=26;


          
%SubPlots=[2,4,6:20,22,24];

minVal = min(min(LV));
maxVal = max(max(LV));
%for iC=1:Nc;
    subplot(sbplt(1),sbplt(2),sbplt(3));
    hold on;
    plot(x,LV(:,Ch(iC)),'color',Color,'LineStyle',LineStyle);
    h=gca;
    set(h,'ylim',[minVal,maxVal]);
    grid on;
    title(ChLab{iC})
    xlabel(xlbl);ylabel(ylbl);
    
%end
%h=gcf;

function  h=topoplot(LV)

%Ch=[1,3,6,8,10,12,14,22,24,26,28,30,40,42,44,46,48,56,58]; %electrodes

chx = [ 53.0000   53.0000  107.0000  121.0000  125.0000  122.0000  114.0000  204.0000  199.0000...
        198.0000  200.0000  206.0000  288.0000  276.0000  271.0000  276.0000  288.0000  341.0000  342.0000].';
    
chy = [ 148.0000  239.0000   71.0000  133.0000  192.0000  255.0000  313.0000   51.0000  119.0000...
        191.0000  259.0000  332.0000   75.0000  128.0000  190.0000  253.0000  305.0000  147.0000  231.0000].';

LV = LV(Ch);
LV=LV(:);

% oldFolder=cd('C:\Users\dionperd\Dropbox\Dionysis\DPtoolbox\Plotting\eegplot');
%oldFolder=cd('C:\Users\CoordAgeEEG\Dropbox\DPtoolbox\Plotting\eegplot');
h=eegplot(double(LV),[chx, chy],1,1,'nearest',[]);
%cd(oldFolder)