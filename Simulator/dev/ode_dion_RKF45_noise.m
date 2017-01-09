function [t, dt, Y]=ode_dion_RKF45_noise(f,TS,et,ndt,idt,xdt,Y0,Q,pin)%,pout0 .... dYdt, d2Ydt2, pout

% Runge-Kutta-Fehlberg Method (RKF45)
% 
% One way to guarantee accuracy in the solution of an I.V.P. is to solve the problem twice using step sizes h and[Graphics:Images/RungeKuttaFehlbergMod_gr_1.gif] and compare answers at the mesh points corresponding to the larger step size.  But this requires a significant amount of computation for the smaller step size and must be repeated if it is determined that the agreement is not good enough. The Runge-Kutta-Fehlberg method (denoted RKF45) is one way to try to resolve this problem.  It has a procedure to determine if the proper step size h is being used.  At each step, two different approximations for the solution are made and compared.  If the two answers are in close agreement, the approximation is accepted. If the two answers do not agree to a specified accuracy, the step size is reduced.  If the answers agree to more significant digits than required, the step size is increased.
% 
% Each Runge-Kutta-Fehlberg step requires the use of the following six values:
%     
% k0 = f(xi, yi)
% 
% k1 = f(xi + 1/4 h, yi + 1/4 k0h)
% 
% k2 = f(xi + 3/8 h, yi + (3/32 k0 + 9/32 k1)h)
% 
% k3 = f(xi + 12/13 h, yi + (1932/2197 k0 – 7200/2197 k1 + 7296/2197 k2)h)
% 
% k4 = f(xi + h, yi + (439/216 k0 - 8 k1 +3680/513 k2 -845/4104 k3)h)
% 
% k5 = f(xi + 1/2 h, yi + (-8/27 k0 + 2 k1 - 3544/2565 k2 + 1859/4104 k3 - 11/40 k4)h)
% 
% 
% Then an approximation to the solution of the I.V.P. is made using a Runge-Kutta method of order 4:
% 
% yi+1 = yi + (25/216 k0 + 1408/2565 k2 + 2197/4104 k3 – 1/5 k4)h
% 
% 
% And a better value for the solution is determined using a Runge-Kutta method of order 5:
% 
% zi+1 = zi + (16/135 k0 + 6656/12825 k2 + 28561/56430 k3 – 9/50 k4 + 2/55 k5)h
% 
% The rest has been Denis-modified:
%--------------------------------------------------------------------------
% The optimal step size sh can be determined by multiplying the scalar s times the current step size h. The scalar s is
% 
% hnew = hold(ehold/(2|z(i+1) – y(i+1)|))^(1/4) = (hold^1.25)* (e/2)^0.25* (|z(i+1) – y(i+1)|)^(-0.25)
% 
% where e is the specified error control tolerance.
%--------------------------------------------------------------------------
%
%New version used:
%  if er<=et increase time step as: dt(it+1) = dt(it) * (e/error)^0.01;
%  if er >et decrease time step as: dt(it+1)   = dt(it) * (e/error)^(1.25/erh(it));
%where "er" is the absolute error per time step, "et" is teh error tolerance, "erh" is the number of
%error hits on this specific iteration, "it" is the iteration index
%For my purposes I choose to recalculate the last time point if the error is bigger than the error tolerance, and to adapt the time step starting from the next time point if the error is too small.

%Gaussian noise of strength Q is added as sqrt(Q)*randn(size(K1))*dt after the solution
%has been selected

%INPUT:
%f: an autonomous vectorized differential equation of dimension D to be integrated
%TS: time span of the integration
%ndt: minimum time step->default ndt0=0.0005
%idt: initial time step->default idt0=0.01
%mdt: maximum time step->default xdt0=0.5
%Y0: initial condition,
%et: error tolerance->default et0=10^(-4);
%Q: noise strength
%p: a cell of parameters of f that might also introduce non-autonomous dependances


%OUTPUT:
%t: time points of the integration
%Y: solution
%dt: the time steps used in the integration
%er: acceptable error of each iteration
%err: arhcive of excessive error hits of each and of the specific iteration
%erh: number of excessive error hits per iteration
%---------------------------------------------------------------------------------------------------


%A noise constants we are going to use:
Q=sqrt(Q);

%Calculate the dimensionalitof the system...
D=length(Y0);
%Dpout=length(pout0);
%...and of its parameters
% Dp=length(p);
%Dimensionality of the statistics:
Ds=1;

%Default values:
ndt0=0.0001;
idt0=0.01;
xdt0=1;
et0=10^(-4);


%Make sure everything is OK. If not give a message and choose the default
%values or exit with an error (if TS error)

%Time span values TS must be monotonously increasing or decreasing
t0=TS(1);                      %Initial time point
tf=TS(end);                    %Final time point
dTS=TS(2:end)-TS(1:end-1);     %Differences of time points
tdir = sign(tf - t0);          %Time direction (positive for increasing times)
if any(tdir*dTS <= 0)
  error('The entries in tspan must strictly increase or decrease!');
end
%Minimum time step has to be smaller than the maximum one
if (ndt>=xdt)
    ndt=ndt0;
    xdt=xdt0;
    disp(['Minimum and maximum time steps changed to the default values ndt0=',num2str(ndt0),' and xdt0=',num2str(xdt0),' because you gave xdt<=ndt values!'])
end
%Initial time step has to be smaller in the [xdt ndt] interval
if (idt<ndt)||(idt>xdt)
    idt=idt0;
    disp(['Initial time step changed to the default value idt0=',num2str(idt0),' because you gave a idt<ndt or idt>xdt value!'])
end
%Error tolerance et has to be positive
if (et<=0)
    et=et0;
    disp(['Error tolerance changed to the default value et0=',num2str(et0),' because you gave a non positive value!'])
end


%A constant for the time step adaptation that we are going to use
% decr=2^(-0.25)*(et^0.25);
% incr=1/2*decr;

%More constants of the Runge-Kutta Algorithm
a  =[1/4 3/32  1932/2197  439/216     -8/27];
b  =[    9/32 -7200/2197   -8          2];
c  =[          7296/2197 3680/513  -3544/2565];
d  =[                    -845/4104  1859/4104];
e  =                                 -11/40;
r4 =[25/216  1408/2565  2197/4104  -1/5];
r5 =[16/135 6656/12825 28561/56430 -9/50 2/55];


%Make an estimation about the total number of iterations, based on the
%initial time step and on the time span
est_it=(tf-t0)/idt;
max_it=10*est_it;
max_er_it=10;
%Allocate memory
Y=zeros(est_it,D);
dYdt=zeros(est_it-1,D);
d2Ydt2=zeros(est_it-2,D);
%pout=zeros(est_it,Dpout);
t=zeros(est_it,1);
dt=zeros(est_it,1);
er=zeros(est_it,1);
erh=zeros(est_it,Ds);
err=zeros(est_it*max_er_it,2);
erC=0;
err=[];

%Initial conditions(5/s(it))*0.25
Y(1,:)=Y0;
%pout(1,:)=pout0;
t(1)=t0;
dt(1)=idt;
it=0;             %actual iterations

disp('-----------------------')
disp('Integration parameters:')
disp(['Error tolerance = ',num2str(et)]);
disp(['Minimum time step =',num2str(ndt)]);
disp(['Maximum time step = ',num2str(xdt)]);
disp(['Initial time step = ',num2str(idt)]);
disp(['Maximum iterations = ',num2str(max_it)]);
disp(['Maximum acceptable errors per iteration =',num2str(max_er_it)]);

error_stop_flag=0;
%Integration
while (t<tf)

    %Let us know the iteration number
    if (mod(it,5000)==0)
        disp(['it = ',num2str(it)])
    end
    if (it>max_it)
        disp(['it = ',num2str(it)])
        %error('Maximum iterations exceeded');
        disp('Maximum iterations exceeded');
        break;
    end
    
    it=it+1;
    repeat=1;

    while (repeat==1)

        if (it>length(erh))
                erh(it)=0;
        end
        if (erh(it)>max_er_it)
%             keyboard
            disp(['t = ',num2str(t(it))])
            disp(['it = ',num2str(it)])
            disp(['dt = ',num2str(dt(it))])
            disp(['Error hits = ',num2str(erh(it))])
            disp(['Error = ',num2str(er(it))])
            %error('System stayed too many times on the same iteration without the error being reduced below error tolerance');
            disp('System stayed too many times on the same iteration without the error being reduced below error tolerance');
            error_stop_flag=1;
            break;
        end

        if (it>2)
            dY=dYdt(it-1,:);
            d2Y=d2Ydt2(it-2,:);
        elseif (it==2)
            dY=dYdt(it-1,:);
            d2Y=zeros(1,D);
        else
            dY=zeros(1,D);
            d2Y=zeros(1,D);
        end


        [k0]  = f(Y(it,:),dY,d2Y,pin);%pout_temp

        [k1] = f(Y(it,:) + a(1)*k0*dt(it)*tdir,dY,d2Y,pin);

        [k2] = f(Y(it,:) + (a(2)*k0 + b(1)*k1)*dt(it)*tdir,dY,d2Y,pin);

        [k3] = f(Y(it,:) + (a(3)*k0 + b(2)*k1 + c(1)*k2)*dt(it)*tdir,dY,d2Y,pin);

        [k4] = f(Y(it,:) + (a(4)*k0 + b(3)*k1 + c(2)*k2 + d(1)*k3)*dt(it)*tdir,dY,d2Y,pin);

        [k5] = f(Y(it,:) + (a(5)*k0 + b(4)*k1 + c(3)*k2 + d(2)*k3 + e*k4)*dt(it)*tdir,dY,d2Y,pin);

        %RK 4th order
        Y4 = Y(it,:) + (r4(1)*k0 + r4(2)*k2 + r4(3)*k3 + r4(4)*k4)*dt(it)*tdir;

        %RK 5th order
        Y5 = Y(it,:) + (r5(1)*k0 + r5(2)*k2 + r5(3)*k3 + r5(4)*k4 + r5(5)*k5)*dt(it)*tdir;
 
        %Calculate error
        er(it)=sqrt(sum((Y5-Y4).^2)); %/sqrt(sum(Y4(it+1,:).^2)), Normalization?...

        if (er(it)<=et)

            %Accept the solution of 4th order and add noise
            Y(it+1,:)=Y4+Q*randn(1,D)*dt(it)*tdir;
            %Extract the output parameters
            %pout(it+1,:) = pout_temp; 
            %Calculate first snd second derivatives of the solution
            dYdt(it,:) = (Y(it+1,:)-Y(it,:))/dt(it);
            if (it>1)
                d2Ydt2(it-1,:) = (dYdt(it,:)-dYdt(it-1,:))/dt(it);
            end
            %Update time step only if there was no excessive error at this
            %time step
            %dt(it+1)= incr * dt(it)^1.25 * er(it)^(-0.25);
  
            if (erh(it)==0)
                dt(it+1)= dt(it) * (et/er(it))^0.01;

                %Keep dt within limits
                if (dt(it+1)>xdt)
                    dt(it+1)=xdt;
                end
            else
                dt(it+1)=dt(it);
            end
            %Update time
            t(it+1)=t(it)+dt(it)*tdir;

            %Change flag in order to advance to the next time step
            repeat=0;

        else

            %Signal this fact for the statistics
            %Record of excessive errors and their timings
            erC=erC+1;
            err(erC,:)=[t(it),er(it)];

            if (it>length(erh))
                erh(it)=0;
            end
            erh(it)=erh(it)+1;

            %Decrease time step in order to repeat the calculation
            %dt(it) = decr * dt(it)^1.25 * er(it)^(-0.25);
            dt(it)= dt(it) * (et/er(it))^(1.25);%/erh(it)

            %Keep dt within limits
            if (dt(it)<ndt)
                %disp(['Minimum time step reached! dt= ', num2str(dt(it)), ', it = ', num2str(it),', erh = ', num2str(erh(it)),', er = ', num2str(er(it))])
                dt(it)=ndt;
            end
        end
    end
    
    if (error_stop_flag==1)
        it=it-1;
        break;
    end
end

%Force all vectors to have the same length:
it=it+1;
er(it)=0;
erh(it)=0;
if (it<est_it)
    t=t(1:it);
    Y=Y(1:it,:);
    %pout=pout(1:it,:);
    dt=dt(1:it);
    er=er(1:it);
    erh=erh(1:it);
end

if (erC<est_it*max_er_it)
    err=err(1:erC,:);
end

% figure(10)
% subplot(4,1,1);
% title('Integration statistics');hold on;
% axis([t0-0.1 tf+0.1 0 max(err(:,2))]); hold on;
% plot(err(:,1),err(:,2),'-ro'); grid on;hold on;
% ylabel('Excessive absolute error'); hold on;
% subplot(4,1,2);
% axis([t0-0.1 tf+0.1 0 max(er)]); hold on;
% plot(t,er,'-bo'); grid on;hold on;%
% ylabel('Acceptable absolute error'); hold on;
% subplot(4,1,3);
% axis([t0-0.1 tf+0.1 0 max(erh)+0.1]); hold on;
% plot(t,erh,'bo'); grid on;hold on;%
% ylabel('Errors per iteration'); hold on;
% subplot(4,1,4);
% axis([t0-0.1 tf+0.1 ndt max(dt)+ndt]); hold on;
% plot(t,dt,'-bo');grid on;hold on;
% xlabel('Time');ylabel('dt'); hold on;

disp('-----------------------')
disp('Integration statistics:')
disp(['Minimum dt = ',num2str(min(dt))]);
disp(['Mean dt = ',num2str(mean(dt))]);
disp(['Maximum dt = ',num2str(max(dt))]);
disp(['Estimated iterations = ',num2str(est_it)]);
iters=length(dt);
disp(['Total iterations = ',num2str(iters)]);
tots=sum(erh);
disp(['Total false calculations = ',num2str(tots)]);
disp(['Total calculations = ',num2str(tots+iters)]);

