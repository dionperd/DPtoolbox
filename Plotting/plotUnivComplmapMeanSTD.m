function h=plotUnivComplmapMeanSTD(M,STD,Ns,Ng,Nc,Groups,Conds,method,measure,cfg)

if strcmpi(method,'STDerr') 
    h = figure('name',[measure, ': standard error'],'units','normalized','outerposition',[0 0 1 1]);
else
    h = figure('name',[measure, ': mean'],'units','normalized','outerposition',[0 0 1 1]);
end

if strcmpi(measure,'MSE')
    x=cfg.MSE.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'MSE';
    for iG = 1:Ng;
        for iC = 1:Nc;
            subplot(Ng,Nc,Nc*(iG-1)+iC)
            plotElectrodes(x,M{iG,iC},STD{iG,iC},Ns,xlbl,ylbl,method);
            hold on;
            title([Groups{iG},' - ',Conds{iC}])
            hold off;
        end
    end
    
elseif strcmpi(measure,'MSEn')
    x=cfg.MSEn.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'MSEn';
    for iG = 1:Ng;
        for iC = 1:Nc;
            subplot(Ng,Nc,Nc*(iG-1)+iC)
            plotElectrodes(x,M{iG,iC},STD{iG,iC},Ns,xlbl,ylbl,method);
            hold on;
            title([Groups{iG},' - ',Conds{iC}])
            hold off;
        end
    end
    
elseif strcmpi(measure,'CMSE')
    x=cfg.CMSE.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'CMSE';
    for iG = 1:Ng;
        for iC = 1:Nc;
            subplot(Ng,Nc,Nc*(iG-1)+iC)
            plotElectrodes(x,M{iG,iC},STD{iG,iC},Ns,xlbl,ylbl,method);
            hold on;
            title([Groups{iG},' - ',Conds{iC}])
            hold off;
        end
    end
    
elseif strcmpi(measure,'CMSEn')
    x=cfg.CMSEn.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'CMSEn';
    for iG = 1:Ng;
        for iC = 1:Nc;
            subplot(Ng,Nc,Nc*(iG-1)+iC)
            plotElectrodes(x,M{iG,iC},STD{iG,iC},Ns(iG),xlbl,ylbl,method);
            hold on;
            title([Groups{iG},' - ',Conds{iC}])
            hold off;
        end
    end
    
elseif strcmpi(measure,'LZ')
    x=cfg.LZ.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'LZ';
    for iG = 1:Ng;
        for iC = 1:Nc;
            subplot(Ng,Nc,Nc*(iG-1)+iC)
            plotElectrodes(x,M{iG,iC},STD{iG,iC},Ns(iG),xlbl,ylbl,method);
            hold on;
            title([Groups{iG},' - ',Conds{iC}])
            hold off;
        end
    end
    
elseif strcmpi(measure,'LZn')
    x=cfg.LZn.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'LZn';
    for iG = 1:Ng;
        for iC = 1:Nc;
            subplot(Ng,Nc,Nc*(iG-1)+iC)
            plotElectrodes(x,M{iG,iC},STD{iG,iC},Ns(iG),xlbl,ylbl,method);
            hold on;
            title([Groups{iG},' - ',Conds{iC}])
            hold off;
        end
    end
    
elseif strcmpi(measure,'STD')
    x=cfg.MSE.scales/cfg.fs;
    xlbl = 'scales (sec)';
    ylbl = 'STD';
    for iG = 1:Ng;
        for iC = 1:Nc;
            subplot(Ng,Nc,Nc*(iG-1)+iC)
            plotElectrodes(x,M{iG,iC},STD{iG,iC},Ns(iG),xlbl,ylbl,method);
            hold on;
            title([Groups{iG},' - ',Conds{iC}])
            hold off;
        end
    end
    
elseif strcmpi(measure,'logF')
    x=log(cfg.DFA.scales/cfg.fs);
    xlbl = 'log(scales) (log(sec))';
    ylbl = 'logF';
    for iG = 1:Ng;
        for iC = 1:Nc;
            subplot(Ng,Nc,Nc*(iG-1)+iC)
            plotElectrodes(x,M{iG,iC},STD{iG,iC},Ns(iG),xlbl,ylbl,method);
            hold on;
            title([Groups{iG},' - ',Conds{iC}])
            hold off;
        end
    end
    
elseif strcmpi(measure,'logV')
    x=log(cfg.VGR.scales/cfg.fs);
    xlbl = 'log(scales) (log(sec))';
    ylbl = 'log(VAR)';
    for iG = 1:Ng;
        for iC = 1:Nc;
            subplot(Ng,Nc,Nc*(iG-1)+iC)
            plotElectrodes(x,M{iG,iC},STD{iG,iC},Ns(iG),xlbl,ylbl,method);
            hold on;
            title([Groups{iG},' - ',Conds{iC}])
            hold off;
        end
    end
    
elseif strcmpi(measure,'logP')
    x=cfg.SP.logf;
    xlbl = 'log(f) (log(Hz))';
    ylbl = 'log(POW)';
    for iG = 1:Ng;
        for iC = 1:Nc;
            subplot(Ng,Nc,Nc*(iG-1)+iC)
            plotElectrodes(x,M{iG,iC},STD{iG,iC},Ns(iG),xlbl,ylbl,method);
            hold on;
            title([Groups{iG},' - ',Conds{iC}])
            hold off;
        end
    end
    
elseif strcmpi(measure,'P')
    x=cfg.SP.f(2:end);
    xlbl = 'f (Hz)';
    ylbl = 'POW';
    for iG = 1:Ng;
        for iC = 1:Nc;
            subplot(Ng,Nc,Nc*(iG-1)+iC)
            plotElectrodes(x,M{iG,iC},STD{iG,iC},Ns(iG),xlbl,ylbl,method);
            hold on;
            title([Groups{iG},' - ',Conds{iC}])
            hold off;
        end
    end
    
elseif strcmpi(measure,'PS')||strcmpi(measure,'VS')||strcmpi(measure,'DOF')||strcmpi(measure,'H')||strcmpi(measure,'mMSE')||strcmpi(measure,'mMSEn')
    for iG = 1:Ng;
        for iC = 1:Nc;
            subplot(Ng,Nc,Nc*(iG-1)+iC)
            topoplot(M{iG,iC},STD{iG,iC},Ns(iG),method);
            hold on;
            title([Groups{iG},' - ',Conds{iC}])
            hold off;
        end
    end
    
end


function h=plotElectrodes(x,M,STD,Ns,xlbl,ylbl,method)


%19 channels, Cz=10:
%Ch=[1,3,6,8,10,12,14,22,24,26,28,30,40,42,44,46,48,56,58]; %electrodes
%ChLab = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8','P7','P3','Pz','P4','P8','O1','O2'};

%58 channels, Cz=26:
Ch=[1:58]; %electrodes
ChLab = { 'Fp1', 'Fpz', 'Fp2', 'AF3', 'AF4', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4',...
          'F6', 'F8', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'T7', 'C5', 'C3',...
          'C1', 'Cz', 'C2', 'C4', 'C6', 'T8', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2',...
          'CP4', 'CP6', 'TP8', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8',...
          'PO7', 'PO3', 'PO1', 'POz', 'PO2', 'PO4', 'PO8', 'O1', 'Oz', 'O2'};
      

hold on;

if isempty(STD) 
    minVal = min(M(:));
    maxVal = max(M(:));
    h=imagesc(x,Ch,M.',[minVal, maxVal]); %x,Ch(),,[minVal, maxVal]

else
    if strcmpi(method,'STDerr')
       y = M./STD * sqrt(Ns);
    else
       y=M;
    end
    y = y.';
    minVal = min(y(:));
    maxVal = max(y(:));
%     alphadata = ones(size(STD));
%     inds = aSTD(STD)<opts.th;
%     minSTD = min(STD(:));
%    alphadata(inds) =  opts.Range(1) + (STD(inds) - minSTD)*( diff(opts.Range)/(opts.th-minSTD) );
    h=imagesc(x,Ch,y);%,'AlphaData',alphadata
    set(gca,'YDir','reverse','clim',[minVal,maxVal]);%  

end
set(gca,'YDir','reverse','ytick',Ch,'yticklabel',ChLab);%  
xlabel(xlbl);ylabel(ylbl);
colorbar;
axis tight
hold off;


function  h=topoplot(M,STD,Ns,method)

% %Ch=[1,3,6,8,10,12,14,22,24,26,28,30,40,42,44,46,48,56,58]; %electrodes
% 
% chx = [ 53.0000   53.0000  107.0000  121.0000  125.0000  122.0000  114.0000  204.0000  199.0000...
%         198.0000  200.0000  206.0000  288.0000  276.0000  271.0000  276.0000  288.0000  341.0000  342.0000].';
%     
% chy = [ 148.0000  239.0000   71.0000  133.0000  192.0000  255.0000  313.0000   51.0000  119.0000...
%         191.0000  259.0000  332.0000   75.0000  128.0000  190.0000  253.0000  305.0000  147.0000  231.0000].';
% M = M(Ch);
% Ch=[1:19]; %electrodes

chx = [ 50.9956,   57.9605,   54.4781,   86.9810,   86.9810,  107.8757,  114.8406,  120.6447,  124.1272,  126.4488, 125.2880,  121.8056,  118.3231,  111.3582,  155.4693,  158.9518,  163.5950,  162.4342,  162.4342,  161.2734, 156.6301,  200.7412,  200.7412,  198.4196,  197.2588,  197.2588,  196.0980,  198.4196,  198.4196,  207.7061, 246.0132,  246.0132,  239.0482,  236.7266,  235.5658,  237.8874,  241.3699,  246.0132,  248.3348,  283.1594, 279.6769,  273.8728,  271.5512,  269.2295,  272.7120,  273.8728,  284.3202, 287.8026,  316.8231,  312.1798, 301.7325,  306.3757,  304.0541, 313.3406,  321.4664,  340.0395,  341.2003,  342.3611].';

chy = [  146.3710  193.3387  237.5968,  135.5323  250.2419   70.5000  100.3065  131.0161  163.5323  190.6290 223.1452  254.7581  283.6613  316.1774   88.5645  123.7903  157.2097  190.6290  224.9516  258.3710 294.5000   50.6290   85.8548  119.2742  153.5968  191.5323  227.6613  261.0806  297.2097  334.2419 52.4355   92.1774  120.1774  152.6935  190.6290  223.1452  255.6613  291.7903  323.4032   75.0161 104.8226  125.5968  158.1129  190.6290  220.4355  252.0484  278.2419  305.3387  103.0161  136.4355 166.2419  187.0161  214.1129  243.9194  272.8226  148.1774  187.9194  229.4677 ]';

if isempty(STD)
    x=M;
else
    if strcmpi(method,'STDerr')
        x = M./STD * sqrt(Ns);
    else
        x=M;
    end
end
x = x(:);
minVal = min(x);
maxVal = max(x);
ColorBar = 1;
h=eegplot(double(x),[chx, chy],1,1,'nearest',[],minVal,maxVal,ColorBar);
hold on;
colormap;
hold off;
