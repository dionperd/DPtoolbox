function [C,clim] = DPconVec2conMatBet(C,cfg)

%Time dimension must be the first one!! Give connectivities as row
%matrices!

[Nt, Nloops] = size(C);

if strcmpi(cfg.freqMode,'CFC')
      
    matC = zeros(Nt,cfg.Nf,cfg.Nf,cfg.Nch,cfg.Nch);
    for iL=1:Nloops;
        matC(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,3),cfg.INDs(iL,4)-21) = C(:,iL);
    end

     if (nargout>1)
        minVal = nan(cfg.Nf,cfg.Nf);
        maxVal = nan(cfg.Nf,cfg.Nf);
        
        for iF1 = 1:cfg.Nf;
            if strcmpi(cfg.chanMode,'directed')
                iF2start = 1;
            else
                iF2start = iF1;
            end
            for iF2 = iF2start:cfg.Nf;
                temp = matC(:,iF1,iF2,:,:);
                temp =temp(:);
                temp = temp(temp~=0);
                minVal(iF1,iF2) = min(temp);
                maxVal(iF1,iF2) = max(temp);
            end
        end
        
        clim = {minVal, maxVal};

     end
    
else
     
    matC = zeros(Nt,cfg.Nf,cfg.Nch,cfg.Nch);
    
    for iL=1:Nloops;
        matC(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,3)-21) = C(:,iL);
    end

     if (nargout>1)
        minVal = nan(cfg.Nf);
        maxVal = nan(cfg.Nf);
        
        for iF = 1:cfg.Nf;
            temp = matC(:,iF,:,:);
            temp =temp(:);
            temp = temp(temp~=0);
            minVal(iF) = min(temp);
            maxVal(iF) = max(temp);    
        end
        
        clim = {minVal, maxVal};

     end
     
end

C=matC;
