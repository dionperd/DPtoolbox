function [x, dxdt, t] = integration(D, x0, tau, t0, iters, dt, s_noise)


%Integration method parameters
atol=1;%10^(-3);
rtol=1;%10^(-3);
Integr_Param={atol rtol};

%Time vector formation
t=[t0:dt:iters*dt]';


%Integration loop that runs in C (mex file)
[x dxdt]=integration_loop_RK45(tau, t, s_noise, Integr_Param, x0');


%A bit of data treatment
x=x.';
dxdt=dxdt.';
minX=min(x);
maxX=max(x);
mindX=min(dxdt);
maxdX=max(dxdt);

%Plotting
figure(1)
for i=1:D;
    subplot(D,1,i);
    axis([-0.1 t(iters)+0.1 minX(i)-0.1 maxX(i)+0.1]);hold on;grid on;
    plot(t,x(:,i),'k');hold on;
    ylabel(['x_',num2str(i)]);
    if (i==D)
        xlabel('t');
    end
    hold off;    
end

figure(2)
for i=1:D;
    subplot(D,1,i);
    axis([-0.1 t(iters)+0.1 mindX(i)-0.1 maxdX(i)+0.1]);hold on;grid on;
    plot(t,dxdt(:,i),'k');hold on;
    ylabel(['dx_',num2str(i),'/dt']);
    if (i==D)
        xlabel('t');
    end
    hold off;    
end

if (D==1)
    figure(3)
    axis([minX-0.1 maxX+0.1 mindX-0.1 maxdX+0.1]);hold on;
    plot(x,dxdt,'k');hold on;
    xlabel('x');ylabel('dx/dt');hold off;
elseif (D==2)
    figure(3)
    axis([minX(1)-0.1 maxX(1)+0.1 minX(2)-0.1 maxX(2)+0.1]);hold on;
    plot(x(:,1),x(:,2),'k');hold on;
    xlabel('x');ylabel('y');hold off;
elseif (D==3)
    figure(3)
    axis([minX(1)-0.1 maxX(1)+0.1 minX(2)-0.1 maxX(2)+0.1 minX(3)-0.1 maxX(3)+0.1]);hold on;
    plot3(x(:,1),x(:,2),x(:,3),'k');hold on;
    xlabel('x');ylabel('y');zlabel('z');hold off;
elseif (D==4)
    figure(3)
    subplot(1,2,1);
    axis([minX(1)-0.1 maxX(1)+0.1 minX(2)-0.1 maxX(2)+0.1 minX(3)-0.1 maxX(3)+0.1]);hold on;
    plot3(x(:,1),x(:,2),x(:,3),'k');hold on;
    xlabel('x_1');ylabel('x_2');zlabel('x_3');hold off;
    subplot(1,2,2);
    axis([minX(3)-0.1 maxX(3)+0.1 minX(4)-0.1 maxX(4)+0.1 minX(1)-0.1 maxX(1)+0.1]);hold on;
    plot3(x(:,3),x(:,4),x(:,1),'k');hold on;
    xlabel('x_3');ylabel('x_4');zlabel('x_1');hold off;
end



