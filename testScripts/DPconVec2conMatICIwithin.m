function [C,clim] = DPconVec2conMatICIwithin(C,cfg)

%Time dimension must be the first one!! Give connectivities as row
%matrices!

[Nt, Nloops] = size(C);
Nloops4 = Nloops/4;
Nloops2 = 2*Nloops4;
Nloops3 = 3*Nloops4;

if strcmpi(cfg.freqMode,'CFC')
    
    matCa = zeros(Nt,cfg.Nf,cfg.Nf,cfg.Nch,cfg.Nch);
    matCb = zeros(Nt,cfg.Nf,cfg.Nf,cfg.Nch,cfg.Nch);
    for iL=1:Nloops4;
        matCa(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,3),cfg.INDs(iL,4)) = C(:,iL);
        matCa(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,4),cfg.INDs(iL,3)) = C(:,iL+Nloops2);
        matCb(:,cfg.INDs(iL+Nloops4,1),cfg.INDs(iL+Nloops4,2),cfg.INDs(iL+Nloops4,3)-21,cfg.INDs(iL+Nloops4,4)-21) = C(:,iL+Nloops4);
        matCb(:,cfg.INDs(iL+Nloops4,1),cfg.INDs(iL++Nloops4,2),cfg.INDs(iL++Nloops4,4)-21,cfg.INDs(iL++Nloops4,3)-21) = C(:,iL+Nloops3);
    end
    
    if (nargout>1)
        minVal = nan(cfg.Nf,cfg.Nf);
        maxVal = nan(cfg.Nf,cfg.Nf);
        
        for iF1 = 1:cfg.Nf;
            for iF2 = 1:cfg.Nf;
                tempA = matCa(:,iF1,iF2,:,:);
                tempB = matCb(:,iF1,iF2,:,:);
                temp = [tempA(:);tempB(:)];
                temp = temp(temp~=0);
                minVal(iF1,iF2) = min(temp);
                maxVal(iF1,iF2) = max(temp);
            end
        end
        
        clim = {minVal, maxVal};

    end
    
else
    
    matCa = zeros(Nt,cfg.Nf,cfg.Nch,cfg.Nch);
    matCb = zeros(Nt,cfg.Nf,cfg.Nch,cfg.Nch);
    for iL=1:Nloops4;
        matCa(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,3)) = C(:,iL);
        matCa(:,cfg.INDs(iL,1),cfg.INDs(iL,3),cfg.INDs(iL,2)) = C(:,iL+Nloops2);
        matCb(:,cfg.INDs(iL+Nloops4,1),cfg.INDs(iL+Nloops4,2)-21,cfg.INDs(iL+Nloops4,3)-21) = C(:,iL+Nloops4);
        matCb(:,cfg.INDs(iL+Nloops4,1),cfg.INDs(iL+Nloops4,3)-21,cfg.INDs(iL+Nloops4,2)-21) = C(:,iL+Nloops3);
    end
    
    if (nargout>1)
        minVal = nan(cfg.Nf);
        maxVal = nan(cfg.Nf);
        
        for iF = 1:cfg.Nf;
            tempA = matCa(:,iF,:,:);
            tempB = matCb(:,iF,:,:);
            temp = [tempA(:);tempB(:)];
            temp = temp(temp~=0);
            minVal(iF) = min(temp);
            maxVal(iF) = max(temp);    
        end
        
        clim = {minVal, maxVal};

    end
       
end

C={matCa, matCb};
