function h=testPh2Ph(method,animation,fc,nm)


cfg.fs=1000;
cfg.fc=fc;
cfg.measures={'IC' 'PLV','PPI','PLI' 'SE', 'CP' , 'MI'};%
cfg.domain='real';
cfg.method = method;
cfg.time=[0:1/cfg.fs:500/cfg.fs]';
cfg.IC.Dphi0=pi/4;
cfg.nm=nm;
cfg.winLen = 0.1;%0.9/cfg.fs;
cfg.cutTails =[0 0];
cfg.timeCalc = cfg.time(1:50:end);
cfg.smoothWinfun = @hanning;
cfg.smoothWinlen = 0.05;
cfg.upsample = 'yes';
% x(:,1,1) = cos(2*pi*cfg.fc(1)*cfg.time);
% x(:,1,2) = x(:,1,1);
% 
% x(:,2,1) = cos(2*pi*cfg.fc(2)*cfg.time - pi/8 + pi/32-pi/16*rand(1001,1));
% x(:,2,2) = cos(2*pi*cfg.fc(2)*cfg.time - pi/4 + pi/8-pi/4*rand(1001,1));
% 
% % x(1:501,2) = cos(2*pi*cfg.fc(2)*cfg.time(1:501) + pi/2 + pi/32-pi/16*rand(501,1));
% x(502:1001,2) = cos(2*pi*cfg.fc(2)*cfg.time(502:1001) - pi/2 + pi/32-pi/16*rand(500,1));

%x(:,2) = cos(2*pi*cfg.fc(2)*cfg.time + 2*pi*rand(1001,1));
% 
%x(:,2) = cos(2*pi*cfg.fc(2)*cfg.time + 2*pi - 4*pi*randn(1001,1));
   

for ii=1:31;  
    InitCondoffset=pi*randn;
    phi1 = 2*pi*cfg.fc(1)*cfg.time + pi/1000*randn(501,1)/cfg.nm(1,1) + InitCondoffset/cfg.nm(1,1);
    phi2(1:501)   = 2*pi*cfg.fc(2)*cfg.time(1:501) + (-pi/6+pi/1000*randn(501,1))/cfg.nm(1,2) + InitCondoffset/cfg.nm(1,2);
    %phi2(1:250)   = 2*pi*cfg.fc(2)*cfg.time(1:250) + (- pi/8 + pi/32-pi/16*rand(250,1))/cfg.nm(1,2) + InitCondoffset/cfg.nm(1,2);
    %phi2(251:501) = 2*pi*cfg.fc(2)*cfg.time(251:501) + (- pi/4 + pi/8-pi/4*rand(251,1))/cfg.nm(1,2) + InitCondoffset/cfg.nm(1,2);
    x(:,1,ii) = cos(phi1);
    x(:,2,ii) = cos(phi2);
end

% for ii=1:101;  
%     x(:,1,ii) = cos(2*pi*cfg.fc(1)*cfg.time + pi/32*randn(501,1)/cfg.nm(2,1));
%     x(1:250,2,ii) = cos(2*pi*cfg.fc(2)*cfg.time(1:250) + (- pi/8 + pi/32-pi/16*rand(250,1))/cfg.nm(2,2));
%     x(251:501,2,ii) = cos(2*pi*cfg.fc(2)*cfg.time(251:501) + (- pi/4 + pi/8-pi/4*rand(251,1))/cfg.nm(2,2));
% end
save(['PC',cfg.method,'.mat'],'x','cfg','method','animation')

[phi cfg method measInds Nmeasures thisfunCommands Ncomnds] = DPphaseCouplingPrepare(x,cfg);
save(['PC',cfg.method,'.mat'],'x','cfg','animation','phi','method', 'measInds', 'Nmeasures', 'thisfunCommands', 'Ncomnds')

[PC, Dphi, cfg, statsRes] = DPphaseCoupling(phi, cfg, method, measInds, Nmeasures, thisfunCommands, Ncomnds);
save(['PC',cfg.method,'.mat'],'x','cfg','animation','phi','method', 'measInds', 'Nmeasures', 'thisfunCommands', 'Ncomnds','PC','Dphi')

h=DPplotPC(x,phi,Dphi,PC,statsRes,cfg,method,measInds, Nmeasures, animation,'test');
