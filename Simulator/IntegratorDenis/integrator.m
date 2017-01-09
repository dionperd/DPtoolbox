function [x dxdt] = integrator(System,param,x0,t,ip)

%Runge Kutta or Euler integrator
%
%Inputs:
%-f: MATLAB handle to the function to be integrated of the form:
%    dxdt=f(t,x,param)
%-param: structure with the parameters of F that may also change in time
%-t: time vector, 
%    i.e. time points where the system is going to be evaluated 
%    (for time step "fix": t(i+1)-t(i)=dt) 
%-x0: state initial condition
%-ip: integraton parameters structure:
%   -ip.mode: "Euler" or "RK_45". Default: "RK_45" (Dormand-Prince method for adaptive 
%       step)
%   -ip.step: "var", "fix". 
%       Default: "fix" ("Euler" integration can be only with "fix" time
%       step)
%   -ip.min_step: minimum step size. Default: 10^(-6)
%   -ip.max_step: maximum step size. Default: 1
%   -ip.s_noise: standard deviation of gaussian noise. Default s_noise=0.
%                (NOTE!: for more accurate results use Runge-Kutta methods
%                for stochastic differential equations!)
%   -ip.err_tol: maximum absolute error allowed per calculation (only for
%                "var" step. Default: 10^(-6).
%
%Outputs:
%-x: solution system's state vector
%-dxdt: function's F evaluations as used for the integration 




%Read or construct integration parameters
if nargin<5
    
    ip=struct('mode',{'RK_45'}, 'step',{'fix'}, 'min_step',10^(-6),...
              'max_step',1, 'err_tol',10^(-6), 's_noise',0);
          
end

if strcmp(ip.step, 'var') 
    ip.step=int32(1);
else
    ip.step=int32(0);
end

if strcmp(ip.mode, 'Euler') 
    ip.mode=int32(0);
else
    ip.mode=int32(1);
end

%Select system to integrate
switch System
    
    case 'excitator'
        System=1;
    
    case 'lorenz'
        System=2;
        
    case 'rossler'
        System=3;
        
    otherwise
        System=0;
            
end

%Integration...

%Make sure all arrays are row arrays
S=size(x0);
if (S(1)<S(2))
    x0=x0.';
end
S=size(t);
if (S(1)<S(2))
    t=t.';
end

%Call the mex/C program to perform the integration
[x dxdt] = integrator_loop(int32(System),param,x0,t,ip);











% %What the MATLAB code does:
% %Preparations...  
% 
% %System's dimension
% D = length(x0);
% 
% %Number of iterations
% N = length(t)-1;
% 
% %Fixed or initial time step
% dt = t(2)-t(1);
%     
% % Preallocate memory
% x = zeros(N+1,D);
% dxdt = zeros(N,D);
% 
% % Assign initial conditions
% x(1,:) = x0;
% 
% 
% 
% 
% 
% if  strcmp(ip.step, 'var')
%     
%     %Butcher's tableau for Dormand-Prince (4,5) method
%     a = [0           0          0           0        0          0;
%         1/5         0          0           0        0          0;
%         3/40        9/40       0           0        0          0;
%         44/45      -56/15      32/9        0        0          0;
%         19372/6561 -25360/2187 64448/6561 -212/729  0          0;
%         9017/3168  -355/33     46732/5247  49/176  -5103/18656 0;
%         35/384	     0	        500/1113	125/192	-2187/6784	11/84];
%     %b5 = [a(7,:) 0] = [35/384 0 500/1113 125/192 -2187/6784 11/84 0]
%     %b4:
%     b = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];
%     c = [0 1/5 3/10 4/5 8/9 1 1];
%     
%     %Initialize error to prevent initial dt from changing (unless outside
%     %limits)
%     error = ip.err_tol;
%     
%     %Integration Loop
%     for ii=1:N;
%         
%         
%         %Set current time point
%         tc = t(ii);
%         
%         %Set current state vector
%         xc = x(ii,:);
%         
%         
%         while (tc<t(ii+1))
%             
%             CALCULATION_DONE = 0;
%             
%             while (CALCULATION_DONE == 0)
%                 
%                 %Adapt step size with the rule dt=dt*(error/tolerance)^(1/5)
%                 %and the constraints dt<t(ii+1)-tc (in order not to overshot the next time point)
%                 %and  min_step<=dt<=max_step
%                 dt = max( [ min([dt*(ip.err_tol/error)^(0.2), t(ii+1)-tc, ip.max_step]), ip.min_step ] );
%                 
%                 %Dormand-Prince method (Runge-Kutta 4-5):
%                 [k1 param] = f(tc,         xc,                                                          param);
%                 [k2 dummy] = f(tc+c(2)*dt, xc+dt*( k1*a(2,1) ),                                         param);
%                 [k3 dummy] = f(tc+c(3)*dt, xc+dt*( k1*a(3,1)+k2*a(3,2) ),                               param);
%                 [k4 dummy] = f(tc+c(4)*dt, xc+dt*( k1*a(4,1)+k2*a(4,2)+k3*a(4,3) ),                     param);
%                 [k5 dummy] = f(tc+c(5)*dt, xc+dt*( k1*a(5,1)+k2*a(5,2)+k3*a(5,3)+k4*a(5,4) ),           param);
%                 [k6 dummy] = f(tc+c(6)*dt, xc+dt*( k1*a(6,1)+k2*a(6,2)+k3*a(6,3)+k4*a(6,4)+k5*a(6,5) ), param);
%                 %[k7 dummy] = f(tc(ii)+c(7)*dt, x(ii,:)+dt*( k1*a(7,1)+k2*a(7,2)+k3*a(7,3)+k4*a(7,4)+k5*a(7,5)+k6*a(7,6) )
%                 %, param);
%                 %but also a(7,2)=0...
%                 a7b5 = k1*a(7,1) +k3*a(7,3)+k4*a(7,4)+k5*a(7,5)+k6*a(7,6);
%                 [k7 dummy] = f(tc+c(7)*dt, xc+dt*( a7b5 )                                             , param);
%                 
%                 %Since, b5=[a(7,:) 0]...
%                 dxdt5 = a7b5;
%                 
%                 %dxdt = b(1)*k1 + b(2)*k2 + b(3)*k3 + b(4)*k4 +b(5)*k5 +b(6)*k6 +b(7)*k7;
%                 %but b(2)=0;
%                 dxdt4 = b(1)*k1  + b(3)*k3 + b(4)*k4 +b(5)*k5 +b(6)*k6 +b(7)*k7;
%                 
%                 %Absolute (euclidean distance) error 
%                 error = 0;
%                 for jj=1:D;
%                     error = error + (dxdt5(jj)-dxdt4(jj))^2;
%                 end
%                 error = sqrt(error)*dt;
%                 
%                 if (error < ip.err_tol)
%                     CALCULATION_DONE = 1;
%                 end
%                 
%                 
%             end
%             
%             
%             %Update current state vector
%             xc = xc + dt*( dxdt4 + ip.s_noise*randn(1,D) );
%             
%             %Advance current time point
%             tc=tc+dt;
%             
%             
%         end
%         
%         %Store next state
%         x(ii+1,:) = xc;
%         
%         %Calculate rate of change
%         dxdt(ii) = ( x(ii+1)-x(ii) ) / ( t(ii+1)-t(ii) );
%         
%     end
%     
%     
%     
%     
% else
%     
%     
%     
%     
%     if strcmp(ip.mode, 'Euler')
%         
%         %Integration Loop
%         for ii=1:N;
%             
%             %Evaluate system's function
%             [dxdt(ii,:) param] = f(t(ii),x(ii,:),param);
%             
%             %Add noise
%             dxdt(ii,:) = dxdt(ii,:)+ip.s_noise*randn(1,D);
%             
%             %Calculate next state
%             x(ii+1,:) = x(ii,:)+dxdt(ii,:)*dt;
%             
%             
%         end
%         
%         
%         
%     else
%         
%         %Butcher's tableau for Runge-Kutta 4rth order method
%         a = [0   0   0;
%             1/2 0   0;
%             0   1/2 0;
%             0   0   1];
%         b = [1/6 1/3 1/3 1/6];
%         c = [0 1/2 1/2 1];
%         
%         %Integration Loop
%         for ii=1:N;
%             
%             %Runge-Kutta 4rth order:
%             [k1 param] = f(t(ii),         x(ii,:),                                    param);
%             [k2 dummy] = f(t(ii)+c(2)*dt, x(ii,:)+dt* k1*a(2,1),                      param);
%             [k3 dummy] = f(t(ii)+c(3)*dt, x(ii,:)+dt*(k1*a(3,1)+k2*a(3,2)),           param);
%             [k4 dummy] = f(t(ii)+c(4)*dt, x(ii,:)+dt*(k1*a(4,1)+k2*a(4,2)+k3*a(4,3)), param);
%             
%             dxdt(ii,:) = b(1)*k1 + b(2)*k2 + b(3)*k3 + b(4)*k4  + ip.s_noise*randn(1,D);
%             
%             %Calculate next state
%             x(ii+1,:) = x(ii,:)+dxdt(ii,:)*dt;
%             
%         end
%         
%         
%         
%     end
%     
% end
% 
% 
% 
% 
% function [f p]=lorenz(t,x,p)
% 
% %For chaos: a=10, b=28, c=8/3
% % p.tau=1;
% % p.sigma=10;
% % p.rho=28;
% % p.beta=8/3;
% %pL=struct('tau',1,'sigma',10,'rho',28,'beta',8/3);
% 
% f(1) = p.sigma*( x(2)-x(1) )/p.tau;
% f(2) = ( x(1)*(p.rho-x(3)) - x(2) )/p.tau;
% f(3) = ( x(1)*x(2) - p.beta*x(3) )/p.tau;
% 
% 
% 
% function [f p]=rossler(t,x,p)
% 
% %For chaos: a=0.1, b=0.1, c=14
% %Rossler studied it for a=0.2, b=0.2, c=5.7
% % p.tau=0.1;
% % p.a=0.1;
% % p.b=0.1;
% % p.c=14;
% %pR=struct('tau',0.1,'a',0.1,'b',0.1,'c',14);
% f(1) = ( -x(2)-x(3) )/p.tau;
% f(2) = ( x(1) + p.a*x(2) )/p.tau;
% f(3) = ( p.b + x(3)*(x(1)-p.c) )/p.tau;
% 
% 
% 
% function [f p]=linear(t,x,p)
% 
% f = -(x-p.x0)/p.tau;

