function [C,clim] = DPconVec2conMatWithin(C,cfg)

%Time dimension must be the first one!! Give connectivities as row
%matrices!

[Nt, Nloops] = size(C);
Nloops2 = Nloops/2;

if strcmpi(cfg.freqMode,'CFC')
      
    matCa = zeros(Nt,cfg.Nf,cfg.Nf,cfg.Nch,cfg.Nch);
    matCb = zeros(Nt,cfg.Nf,cfg.Nf,cfg.Nch,cfg.Nch);
    for iL=1:Nloops2;
        matCa(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,3),cfg.INDs(iL,4)) = C(:,iL);
        matCb(:,cfg.INDs(Nloops2+iL,1),cfg.INDs(Nloops2+iL,2),cfg.INDs(Nloops2+iL,3)-21,cfg.INDs(Nloops2+iL,4)-21) = C(:,Nloops2+iL);
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
    for iL=1:Nloops2;
        matCa(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,3)) = C(:,iL);
        matCb(:,cfg.INDs(Nloops2+iL,1),cfg.INDs(Nloops2+iL,2)-21,cfg.INDs(Nloops2+iL,3)-21) = C(:,Nloops2+iL);
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