function [C,clim] = DPconVec2conMatICIbet(C,cfg)

%Time dimension must be the first one!! Give connectivities as row
%matrices!

[Nt, Nloops] = size(C);
Nloops2 = Nloops/2;

if strcmpi(cfg.freqMode,'CFC')
      
    matC12 = zeros(Nt,cfg.Nf,cfg.Nf,cfg.Nch,cfg.Nch);
    matC21 = zeros(Nt,cfg.Nf,cfg.Nf,cfg.Nch,cfg.Nch);
    for iL=1:Nloops2;
        matC12(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,3),cfg.INDs(iL,4)-21) = C(:,iL);
        matC21(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,3),cfg.INDs(iL,4)-21) = C(:,iL+Nloops2);
    end

   if (nargout>1)
        minVal = nan(cfg.Nf,cfg.Nf);
        maxVal = nan(cfg.Nf,cfg.Nf);
        
        for iF1 = 1:cfg.Nf;
            for iF2 = 1:cfg.Nf;
                temp12 = matC12(:,iF1,iF2,:,:);
                temp21 = matC21(:,iF1,iF2,:,:);
                temp = [temp12(:);temp21(:)];
                temp = temp(temp~=0);
                minVal(iF1,iF2) = min(temp);
                maxVal(iF1,iF2) = max(temp);
            end
        end
        
        clim = {minVal, maxVal};

    end
        
else
     
    matC12 = zeros(Nt,cfg.Nf,cfg.Nch,cfg.Nch);
    matC21 = zeros(Nt,cfg.Nf,cfg.Nch,cfg.Nch);
    for iL=1:Nloops2;
        matC12(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,3)-21) = C(:,iL);
        matC21(:,cfg.INDs(iL,1),cfg.INDs(iL,2),cfg.INDs(iL,3)-21) = C(:,iL+Nloops2);
    end
    
    if (nargout>1)
        minVal = nan(cfg.Nf);
        maxVal = nan(cfg.Nf);
        
        for iF = 1:cfg.Nf;
            temp12 = matC12(:,iF,:,:);
            temp21 = matC21(:,iF,:,:);
            temp = [temp12(:);temp21(:)];
            temp = temp(temp~=0);
            minVal(iF) = min(temp);
            maxVal(iF) = max(temp);    
        end
        
        clim = {minVal, maxVal};

    end

end

C = {matC12, matC21};