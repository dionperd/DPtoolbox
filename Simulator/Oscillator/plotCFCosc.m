function plotCFCosc(t,x)


figure;
subplot(3,1,1)

plot(t,x(:,1).*sin(x(:,2)),'b');hold on;
plot(t,x(:,3).*sin(x(:,4)),'g');hold off;


subplot(3,1,2)

plot(t,x(:,1),'b');hold on;
plot(t,x(:,3),'g');hold off;


subplot(3,1,3)

plot(t,x(:,2),'b');hold on;
plot(t,x(:,4),'g');
plot(t,x(:,4)-x(:,2),'r');hold off;