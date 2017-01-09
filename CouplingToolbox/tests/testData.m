function testData(type,example,fs,fc,time,M)




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
   

for ii=1:101;  
    InitCondoffset=pi*randn;
    phi1 = 2*pi*cfg.fc(1)*cfg.time + pi/32*randn(501,1)/cfg.nm(1,1) + InitCondoffset/cfg.nm(1,1);
    phi2(1:250)   = 2*pi*cfg.fc(2)*cfg.time(1:250) + (- pi/8 + pi/32-pi/16*rand(250,1))/cfg.nm(1,2) + InitCondoffset/cfg.nm(1,2);
    phi2(251:501) = 2*pi*cfg.fc(2)*cfg.time(251:501) + (- pi/4 + pi/8-pi/4*rand(251,1))/cfg.nm(1,2) + InitCondoffset/cfg.nm(1,2);
    x(:,1,ii) = cos(phi1);
    x(:,2,ii) = cos(phi2);
end

