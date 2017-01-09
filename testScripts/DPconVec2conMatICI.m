function [C,clim] = DPconVec2conMatICI(C,cfg)

%Time dimension must be the first one!! Give connectivities as row
%matrices!

[Nt, Nloops] = size(C);

if strcmpi(cfg.freqMode,'CFC')
    
    matC = zeros(Nt,cfg.Nf,cfg.Nf,cfg.Nch,cfg.Nch);
    for iL=1:Nloops/2;
        matC(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,3),cfg.INDs(iL,4)) = C(:,iL);
        matC(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,4),cfg.INDs(iL,3)) = C(:,iL+Nloops/2);
    end
    
    if (nargout>1)
        minVal = nan(cfg.Nf,cfg.Nf);
        maxVal = nan(cfg.Nf,cfg.Nf);
        
        for iF1 = 1:cfg.Nf;
            for iF2 = 1:cfg.Nf;
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
    
    for iL=1:Nloops/2;
        matC(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,3)) = C(:,iL);
        matC(:,cfg.INDs(iL,1),cfg.INDs(iL,3),cfg.INDs(iL,2)) = C(:,iL+Nloops/2);
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
