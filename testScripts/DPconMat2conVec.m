function [C,cfg] = DPconMat2conVec(C,cfg)

[Labels, INDs, fc, nm] = DPconMat2conVectCFG(cfg.Nf,cfg.Nch,cfg.f,cfg.hdr.label,cfg.freqMode,cfg.chanMode);
clear fc nm;

Nloops = size(INDs,1);

%if there is a time dimension as well, put it in the end, assuming it is
%the first one:
Csize = size(C);
NcSize = length(Csize);

if strcmpi(cfg.freqMode,'CFC') || strcmpi(cfg.freqMode,'am2ph')
    
    if NcSize>4 %if there is time dimension (always the first one)
        
        vecC = zeros(Csize(1),Nloops);
        for iL = 1:Nloops;
            vecC(:,iL) = C(:,INDs(iL,1),INDs(iL,2),INDs(iL,3),INDs(iL,4));
        end
    else
        vecC = zeros(1,Nloops);
        for iL = 1:Nloops;
            vecC(1,iL) = C(INDs(iL,1),INDs(iL,2),INDs(iL,3),INDs(iL,4));
        end
    end
    
else
    if NcSize>3 %if there is time dimension (always the first one)
        
        vecC = zeros(Csize(1),Nloops);
        for iL = 1:Nloops;
            vecC(:,iL) = C(:,INDs(iL,1),INDs(iL,2),INDs(iL,3));
        end
    else
        
        vecC = zeros(1,Nloops);
        for iL = 1:Nloops;
            vecC(1,iL) = C(INDs(iL,1),INDs(iL,2),INDs(iL,3));
        end
    end
    
end
C=vecC;
cfg.Labels = Labels;
cfg.INDs = INDs;