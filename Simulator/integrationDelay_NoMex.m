 function [ t, x, pout, dxdt, dx, cfg] = integrationDelay_NoMex(System,p1,t1,x01,dxdt01,cfg1,plotting)
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!WRONG: noise is added for RK, as if it wer Euler!!!!!!!!!!!!!!
%General purpose integrator, "Euler" or "Runge-Kutta 4(5)" for fixed or adaptive time step, with or without noise, for autonomous or non autonomous systems 

global p t x0 dxdt0 cfg;
p=p1;
t=t1;
x0=x01;
dxdt0=dxdt01;
cfg = cfg1;
clear p1 t1 x01 dxdt01 cfg1;


%Input:

%-System: a function handle to the function that is going to be integrated,
%it has to be of the form: [dxdt, pout] = System(t,x,p)

%-p: a structure containing all parameters of the system to be integrated,
%    or function handles for their calculations of the general form @(t,x,dxdt)

%-t: time vector, a column vector of monotonously increasing or decreasing
%real numbers

%-x0: the initial condition, a matrix of real numbers of a dimensions x trials size 
%     it is also used as a history function for t-lag<t(1)
%-dxdt0: the initial condition for the derivative of the system to be used
%        as a history function for t-lag<t(1)
%        a matrix of real numbers of a dimensions x trials size
%        default = [];
%-cfg: a structure with the parameters of the integration:
%     -lags: a vector of lags in time, 
%            monotonously increasing positive real numbers 
%            default = []
%     -method: it has to be one of the following strings: "Euler" or "RK4"
%     (for Runge-Kutta), default is "RK4"
%     -systemType: it has to be one of the following strings: "auto" or
%                 "nonauto", default is "auto"
%     -timestep: it has to be one of the following strings: "fix" or
%                "adapt", default is "fix"
%     -s_noise: the standard deviation of the gaussian random noise, 
%               either per dimension of x, or one for all of them
%     -noise_method: one of the following strings: "add" or
%                    "multi", default is "add"
%     -error_check: one of the following strings: "max" or
%                   "norm", default is "norm"
%     -atol: absolute error, real positive scalar or vector of the system's dimensionality, default is 10^-3
%     -rtol: relative error, real positive scalar or vector of the system's dimensionality, default is 10^-3
%     -ndt: minimum time step, real scalar
%     -mdt: maximum time step, real scalar
%     -idt: fixed or initial time step, real scalar
%     -dt_min_ratio: minimum ratio of PI timestep control, positive real scalar, default is 1/5
%     -dt_max_ratio: maximum ratio of PI timestep control, positive real scalar, default is 5
%     -dt_a: exponent of current error, positive real scalar, default is a=0.7/k, where k is the order of the method (5 for Runge-Kutta-Fehlberg)
%     -dt_b: exponent of previous iteration error, positive real scalar, default is b=0.4/k, where k is the order of the method (5 for Runge-Kutta-Fehlberg)
%     -dt_S: safety parameter for the timestep control, positive real scalr in the interval (0,1), default is 0.9
%     -max_iters: maximum iterations per time step, positive integer scalar, default  is 100
%     -max_tot_steps: maximum total steps, positive integer scalar, default  is 10^6


%Output:

%-t: time vector

%-x: the integrated output, a real matrix of time x dimension x trials, or time x dimension (if there is only one trial) structure

%-dxdt: the first derivative before noise application, a real matrix of time x dimension x trials, or time x dimension (if there is only one trial) structure

%-dx: the increments after noise application, a real matrix of time x dimension x trials, or time x dimension (if there is only one trial) structure

%-pout: the structure of the parameters of the system, in a time x trials
%    format if the system is nonautonomous and the trials are more than 1

%-cfg: the structure with the parameters of the integration extended with
%      the statistics of the integration
%      dt: the time steps used in the integration
%      er: acceptable error of each iteration
%      ers: archive of excessive error hits of each and of the specific iteration
%      erh: number of excessive error hits per iteration


%For tha adaptive time step Runge Kutta method used here:
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
%  if er >et decrease time step as: dt(it)   = dt(it) * (e/error)^(1.25/erh(it));
%where "er" is is the absolute error per time step, "et" is teh error tolerance, "erh" is the number of
%error hits on this specific iteration, "it" is the iteration index
%For my purposes I choose to recalculate the last time point if the error
%is bigger than the error tolerance, and to adapt the time step starting from the next time point if the error is too small.



%Input data validation
funcName='integration_NoMex';

varName='x0';
testX0 = {@(x0)isnumeric(x0),...
          @(x0)isreal(x0),...
          @(x0)length(size(x0))==2 };
param={{},{},{}};
mode=['e','e','e'];
execfun={{},{},{}};
default='';
[~, x0] = DPvalidateData(x0,testX0,param,mode,execfun,default,varName,funcName);
%The number of dimensions and trials
[cfg.D, cfg.Ntr]=size(x0);

if ~isempty(dxdt0)
    varName='dxdt0';
    testDXDT0 = {@(dxdt0)isnumeric(dxdt0),...
        @(dxdt0)isreal(dxdt0),...
        @(dxdt0,x0)(size(dxdt0))==size(x0) };
    param={{},{},{x0}};
    mode=['e','e','e'];
    execfun={{},{},{}};
    default=[];
    [~, dxdt0] = DPvalidateData(dxdt0,testDXDT0,param,mode,execfun,default,varName,funcName);
end

varName='t';
testT = {@(t)isnumeric(t),...
         @(t)isvector(t),...
         @(t)isreal(t),...
         @(t)all(diff(t)>0)||all(diff(t)<0),...
         @(t)size(t,1)>=size(t,2) };
param={{},{},{},{},{}};
mode=['e','e','e','e','w'];
execfun={{},{},{},{},@(t)t'};
default='';
[~, t] = DPvalidateData(t,testT,param,mode,execfun,default,varName,funcName);
%The number of time points
cfg.Nt = length(t);
% %The direction of the integration in time 
difft = diff(t);
tdir = sign(difft(1)); 

varName='System';
testSYSTEM = {@(System)isa(System,'function_handle') };
param={{}};
mode=['e'];
execfun={{}};
default='';
[~, System] = DPvalidateData(System,testSYSTEM,param,mode,execfun,default,varName,funcName);    

varName='p';
testP = {@(p)isstruct(p)||isempty(p) };
param={{}};
mode=['e'];
execfun={{}};
default='';
[~, p] = DPvalidateData(p,testP,param,mode,execfun,default,varName,funcName);   


varName='cfg';
testCFG = { @(cfg)isstruct(cfg) };
param={{}};
mode=['e'];
execfun={{}};
default=createCFG;
[RESULTcfg cfg] = DPvalidateData(cfg,testCFG,param,mode,execfun,default,varName,funcName);   


if RESULTcfg
    
    if isfield(cfg,'lags')
        if ~isempty(cfg.lags)
            varName='lags';
            testLAGS = {@(lags)isnumeric(lags),...
                        @(lags)isvector(lags),...
                        @(lags)isreal(lags),...
                        @(lags)all( diff(lags)>0 )  };
            param={{},{},{},{}};
            mode=['e','e','e','w'];
            execfun={{},{},{},@(lags)unique(lags)};
            default='';
            [~, cfg.lags] = DPvalidateData(cfg.lags,testLAGS,param,mode,execfun,default,varName,funcName);
        end
        
        if (any(cfg.lags==0))
            warning('Getting rid of 0 lag found')
            cfg.lags(cfg.lags==0)=[];
        end
    else
        warning('No vector for lags is given. Setting default cfg.lags = []')
        cfg.lags = [];
    end
    cfg.Nlags = length(cfg.lags);
    
    if isfield(cfg,'method')
        varName='method';
        testMETHOD = {@(method)ischar(method),...
                      @(method)isvector(method),...
                      @(method)strcmpi(method,'RK4')||strcmpi(method,'Euler') };
        param={{},{},{}};
        mode=['e','e','e'];
        execfun={{},{},{}};
        default='RK4';
        [~, cfg.method] = DPvalidateData(cfg.method,testMETHOD,param,mode,execfun,default,varName,funcName);
    else
        cfg.method='RK4';
    end
    
    if isfield(cfg,'systemType')
        varName='systemType';
        testSYSTEMTYPE = {@(systemType)ischar(systemType),...
                          @(systemType)isvector(systemType),...
                          @(systemType)strcmpi(systemType,'auto')||strcmpi(systemType,'nonauto') };
        param={{},{},{}};
        mode=['e','e','e'];
        execfun={{},{},{}};
        default='auto';
        [~, cfg.systemType] = DPvalidateData(cfg.systemType,testSYSTEMTYPE,param,mode,execfun,default,varName,funcName);
    else
        cfg.systemType='auto';
    end

    if isfield(cfg,'s_noise')
        varName='s_noise';
        testSNOISE = {@(s_noise)isnumeric(s_noise),...
                      @(s_noise,D)isscalar(s_noise)|| ( isvector(s_noise)&&(length(s_noise)==D) ),...
                      @(s_noise)isreal(s_noise),...
                      @(s_noise)all(s_noise>=0) };
        param={{},{cfg.D},{},{}};
        mode=['e','e','e','e'];
        execfun={{},{},{},{}};
        default=0;
        [~, cfg.s_noise] = DPvalidateData(cfg.s_noise,testSNOISE,param,mode,execfun,default,varName,funcName);
    else
        cfg.s_noise=0;
    end
    if any(cfg.s_noise>0)
        NOISE=1;
    else
        NOISE=0;
    end
    
    if NOISE
    
        if isfield(cfg,'noise_method')
            varName='noise_method';
            testNOISEMETHOD = {@(noise_method)ischar(noise_method),...
                               @(noise_method)isvector(noise_method),...
                               @(noise_method)strcmpi(noise_method,'add')||strcmpi(noise_method,'multi') };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default='add';
            [~, cfg.noise_method] = DPvalidateData(cfg.noise_method,testNOISEMETHOD,param,mode,execfun,default,varName,funcName);
        else
            cfg.noise_method='add';
        end
        
    end
    
    
    if isfield(cfg,'timestep')
        varName='timestep';
        testTIMESTEP = {@(timestep)ischar(timestep),...
                        @(timestep)isvector(timestep),...
                        @(timestep)strcmpi(timestep,'fix')||strcmpi(timestep,'adapt') };
        param={{},{},{}};
        mode=['e','e','e'];
        execfun={{},{},{}};
        default='fix';
        [~, cfg.timestep] = DPvalidateData(cfg.timestep,testTIMESTEP,param,mode,execfun,default,varName,funcName);
    else
        cfg.timestep='fix';
    end
    
    if strcmpi(cfg.timestep,'adapt')
        
        
        order=5; %For Runge-Kutta-Fehlberg (only method implemented so far) the order of the method is 5
        
        if isfield(cfg,'error_check')
            varName='error_check';
            testERRORCHECK = {@(error_check)ischar(error_check),...
                              @(error_check)isvector(error_check),...
                              @(error_check)strcmpi(error_check,'norm')||strcmpi(error_check,'max') };
            param={{},{},{}};
            mode=['e','e','e'];
            execfun={{},{},{}};
            default='norm';
            [~, cfg.error_check] = DPvalidateData(cfg.error_check,testERRORCHECK,param,mode,execfun,default,varName,funcName);
        else
            cfg.error_check='norm';
        end
        
        if isfield(cfg,'atol')
            varName='atol';
            testATOL = {@(atol)isnumeric(atol),...
                        @(atol)isscalar(atol),...
                        @(atol)isreal(atol),...
                        @(atol)atol>0 };
            param={{},{},{},{}};
            mode=['e','w','e','w'];
            execfun={{},@(atol)atol(1),{},@(atol)-atol};
            default=10^-3;
            [~, cfg.atol] = DPvalidateData(cfg.atol,testATOL,param,mode,execfun,default,varName,funcName);
        else
            cfg.atol=10^-3;
        end
        
        if isfield(cfg,'rtol')
            varName='rtol';
            testRTOL = {@(rtol)isnumeric(rtol),...
                        @(rtol)isscalar(rtol),...
                        @(rtol)isreal(rtol),...
                        @(rtol)rtol>0 };
            param={{},{},{},{}};
            mode=['e','w','e','w'];
            execfun={{},@(rtol)rtol(1),{},@(rtol)-rtol};
            default=10^-3;
            [~, cfg.rtol] = DPvalidateData(cfg.rtol,testRTOL,param,mode,execfun,default,varName,funcName);
        else
            cfg.atol=10^-3;
        end
        
        %Default values for timesteps
        difft=abs(difft);
        dt_min=min(difft)*tdir;
        dt_max=max(difft)*tdir;
        dt_in = (dt_max + dt_min)/2;
        
        if isfield(cfg,'idt')
            varName='idt';
            testIDT = {@(idt)isnumeric(idt),...
                       @(idt)isscalar(idt),...
                       @(idt)isreal(idt),...
                       @(idt,tdir)sign(idt)==tdir };
            param={{},{},{},{tdir}};
            mode=['e','w','e','w'];
            execfun={{},@(idt)idt(1),{},@(idt)-idt};
            default=dt_in;
            [~, cfg.idt] = DPvalidateData(cfg.idt,testIDT,param,mode,execfun,default,varName,funcName);
        else
            cfg.idt=dt_in;
        end
        
        if isfield(cfg,'ndt')
            varName='ndt';
            testNDT = {@(ndt)isnumeric(ndt),...
                       @(ndt)isscalar(ndt),...
                       @(ndt)isreal(ndt),...
                       @(ndt,tdir)sign(ndt)==tdir,...
                       @(ndt,idt)abs(ndt)<abs(idt)};
            param={{},{},{},{tdir},{cfg.idt}};
            mode=['e','w','e','w','e'];
            execfun={{},@(ndt)ndt(1),{},@(ndt)-ndt,{}};
            default=dt_min;
            [~, cfg.ndt] = DPvalidateData(cfg.ndt,testNDT,param,mode,execfun,default,varName,funcName);
        else
            cfg.ndt=dt_min;
        end
        
        if isfield(cfg,'mdt')
            varName='mdt';
            testMDT = {@(mdt)isnumeric(mdt),...
                       @(mdt)isscalar(mdt),...
                       @(mdt)isreal(mdt),...
                       @(mdt,tdir)sign(mdt)==tdir,...
                       @(ndt,idt)abs(mdt)>abs(idt)};
            param={{},{},{},{tdir},{cfg.idt}};
            mode=['e','w','e','w','e'];
            execfun={{},@(mdt)mdt(1),{},@(mdt)-mdt,{}};
            default=dt_max;
            [~, cfg.mdt] = DPvalidateData(cfg.mdt,testMDT,param,mode,execfun,default,varName,funcName);
        else
            cfg.idt=dt_max;
        end
        
        if isfield(cfg,'dt_min_ratio')
            varName='dt_min_ratio';
            testDTMINRATIO = {@(dt_min_ratio)isnumeric(dt_min_ratio),...
                              @(dt_min_ratio)isscalar(dt_min_ratio),...
                              @(dt_min_ratio)isreal(dt_min_ratio),...
                              @(dt_min_ratio)(dt_min_ratio>0)&&(dt_min_ratio<1) };
            param={{},{},{},{}};
            mode=['e','w','e','e'];
            execfun={{},@(dt_min_ratio)dt_min_ratio(1),{},{}};
            default=0.2;
            [~, cfg.dt_min_ratio] = DPvalidateData(cfg.dt_min_ratio,testDTMINRATIO,param,mode,execfun,default,varName,funcName);
        else
            cfg.dt_min_ratio=0.2;
        end
        
        if isfield(cfg,'dt_max_ratio')
            varName='dt_max_ratio';
            testDTMAXRATIO = {@(dt_max_ratio)isnumeric(dt_max_ratio),...
                              @(dt_max_ratio)isscalar(dt_max_ratio),...
                              @(dt_max_ratio)isreal(dt_max_ratio),...
                              @(dt_max_ratio)dt_max_ratio>1 };
            param={{},{},{},{}};
            mode=['e','w','e','e'];
            execfun={{},@(dt_max_ratio)dt_max_ratio(1),{},{}};
            default=5;
            [~, cfg.dt_max_ratio] = DPvalidateData(cfg.dt_min_ratio,testDTMAXRATIO,param,mode,execfun,default,varName,funcName);
        else
            cfg.dt_max_ratio=5;
        end
        
        if isfield(cfg,'dt_a')
            varName='dt_a';
            testDTA = {@(dt_a)isnumeric(dt_a),...
                       @(dt_a)isscalar(dt_a),...
                       @(dt_a)isreal(dt_a),...
                       @(dt_a)(dt_a<0)&&(dt_a>=-1) };
            param={{},{},{},{}};
            mode=['e','w','e','e'];
            execfun={{},@(dt_a)dt_a(1),{},{}};
            default=-0.7/order;
            [~, cfg.dt_a] = DPvalidateData(cfg.dt_a,testDTA,param,mode,execfun,default,varName,funcName);
        else
            cfg.dt_a=-0.7/order;
        end
        
        if isfield(cfg,'dt_b')
            varName='dt_b';
            testDTB = {@(dt_b)isnumeric(dt_b),...
                       @(dt_b)isscalar(dt_b),...
                       @(dt_b)isreal(dt_b),...
                       @(dt_b)(dt_b>0)&&(dt_b<=1) };
            param={{},{},{},{}};
            mode=['e','w','e','e'];
            execfun={{},@(dt_b)dt_b(1),{},{}};
            default=0.4/order;
            [~, cfg.dt_b] = DPvalidateData(cfg.dt_b,testDTB,param,mode,execfun,default,varName,funcName);
        else
            cfg.dt_b=0.4/order;
        end
        
        if isfield(cfg,'dt_S')
            varName='dt_S';
            testDTS = {@(dt_S)isnumeric(dt_S),...
                       @(dt_S)isscalar(dt_S),...
                       @(dt_S)isreal(dt_S),...
                       @(dt_S)(dt_S>0)&&(dt_S<=1) };
            param={{},{},{},{}};
            mode=['e','w','e','e'];
            execfun={{},@(dt_S)dt_S(1),{},{}};
            default=0.9;
            [~, cfg.dt_S] = DPvalidateData(cfg.dt_S,testDTS,param,mode,execfun,default,varName,funcName);
        else
            cfg.dt_S=0.9;
        end
        
        if isfield(cfg,'max_iters')
            varName='max_iters';
            testMAXITERS = {@(max_iters)isnumeric(max_iters),...
                            @(max_iters)isscalar(max_iters),...
                            @(max_iters)isreal(max_iters),...
                            @(max_iters)round(max_iters)==max_iters,...,...
                            @(max_iters)max_iters>0 };
            param={{},{},{},{},{}};
            mode=['e','w','e','w','w'];
            execfun={{},@(max_iters)max_iters(1),{},@(max_iters)round(max_iters),@(max_iters)-max_iters};
            default=10;
            [~, cfg.max_iters] = DPvalidateData(cfg.max_iters,testMAXITERS,param,mode,execfun,default,varName,funcName);
        else
            cfg.max_iters=10;
        end
        
        if isfield(cfg,'max_tot_steps')
            varName='max_tot_steps';
            testMAXTOTSTEPS = {@(max_tot_steps)isnumeric(max_tot_steps),...
                               @(max_tot_steps)isscalar(max_tot_steps),...
                               @(max_tot_steps)isreal(max_tot_steps),...
                               @(max_tot_steps)round(max_tot_steps)==max_tot_steps,...,...
                               @(max_tot_steps)max_tot_steps>0 };
            param={{},{},{},{},{}};
            mode=['e','w','e','w','w'];
            execfun={{},@(max_tot_steps)max_tot_steps(1),{},@(max_tot_steps)round(max_tot_steps),@(max_tot_steps)-max_tot_steps};
            default=10^6;
            [~, cfg.max_tot_steps] = DPvalidateData(cfg.max_tot_steps,testMAXTOTSTEPS,param,mode,execfun,default,varName,funcName);
        else
            cfg.max_tot_steps=10^6;
        end
        
    end


end




%Initialization
global pout x dxdt;
pout=[];
x=zeros(cfg.Nt,cfg.D,cfg.Ntr);
dxdt=zeros(cfg.Nt,cfg.D,cfg.Ntr);


if NOISE
    try
        RandStream.setGlobalStream...
            (RandStream('mt19937ar','seed',sum(100*clock)));
    catch
        RandStream.setDefaultStream ...
            (RandStream('mt19937ar','seed',sum(100*clock)));
    end
    
    dx=zeros(cfg.Nt,cfg.D,cfg.Ntr);
    switch cfg.noise_method
        case 'multi'
            noiseFun=@multi_noise;
        otherwise
            noiseFun=@add_noise;    
    end
end


if strcmpi(cfg.timestep,'fix')
    
    %Permute the resulting vectors so that time is the last index
    x=permute(x,[2,3,1]);
    dxdt=permute(dxdt,[2,3,1]);
    if NOISE
        dx=permute(dx,[2,3,1]);
        %...and make s_noise to have the same size as x of a given time point
        %...i.e. [D,Ntr], unless it is a scalar
        if numel(cfg.s_noise)~=1
            if size(cfg.s_noise,1) == cfg.D
                cfg.s_noise = repmat(cfg.s_noise,1,cfg.Ntr);
            else
                cfg.s_noise = repmat(cfg.s_noise.',1,cfg.Ntr);
            end
        end
    end
    
    %Choose integrator
    if strcmpi(cfg.method,'RK4')
        Integrator = @RK4;
    else
        Integrator = @Euler;
    end

    %Time steps are fixed and given by the difference vector of time
    cfg.dt = diff(t);
    
    %Initialization
    x(:,:,1)=x0;
  
    
    %In this case we calculate all trials at the same time
    
    %Integration loop
    for ii=2:cfg.Nt;
        
        if (cfg.Nlags>0)
            [ tlag,xlag,dxdtlag,plag ] = set_lagsFIX(ii);
        else %...if there are no lags...
            xlag = [];
            dxdtlag = [];
            tlag = [];
            plag = [];
        end
        
        %Calculate first derivative
        if strcmpi(cfg.systemType,'auto')
            [dxdt(:,:,ii-1), ~] = Integrator(t(ii-1),squeeze(x(:,:,ii-1)),tlag,xlag,dxdtlag,cfg.dt(ii-1),System,p,plag);
        else
            [dxdt(:,:,ii-1), pout(ii-1,:)] = Integrator(t(ii-1),squeeze(x(:,:,ii-1)),tlag,xlag,dxdtlag,cfg.dt(ii-1),System,p,plag);
        end
        
        cdx= cfg.dt(ii-1)*dxdt(:,:,ii-1);
        if NOISE
            dx(:,:,ii-1) = cdx + ...
                noiseFun(t(ii-1),dxdt(:,:,ii-1),cfg.s_noise,cfg.dt(ii-1));
            cdx = dx(:,:,ii-1);
        end
        x(:,:,ii) = x(:,:,ii-1) + cdx;
    end
    
    if (cfg.Nlags>0)
            [ tlag,xlag,dxdtlag,plag ] = set_lagsFIX(ii+1);
        else %...if there are no lags...
            xlag = [];
            dxdtlag = [];
            tlag = [];
            plag = [];
    end
        
    %Calculate first derivatives for the last time point
    if strcmpi(cfg.systemType,'auto')
        [dxdt(:,:,ii), ~] = Integrator(t(ii),squeeze(x(:,:,ii)),tlag,xlag,dxdtlag,cfg.dt(ii-1),System,p,plag);
    else
        [dxdt(:,:,ii), pout(ii,:)] = Integrator(t(ii),squeeze(x(:,:,ii)),tlag,xlag,dxdtlag,cfg.dt(ii-1),System,p,plag);
    end
    if NOISE
        dx(:,:,ii) = dxdt(:,:,ii)*cfg.dt(ii-1) + noiseFun(t(ii),dxdt(:,:,ii),cfg.s_noise,cfg.dt(ii-1));
    end
    
    %Permute the resulting vectors so that time is the first index
    x=permute(x,[3,1,2]);
    dxdt=permute(dxdt,[3,1,2]);
    if NOISE
        dx=permute(dx,[3,1,2]);
    end
    
else
    
    %Choose integrator
    Integrator = @RK45;
    
    %Make s_noise a row vector if it is not a scalar
    if NOISE
        if numel(cfg.s_noise)~=1
            if size(cfg.s_noise,2) ~= cfg.D
                cfg.s_noise = cfg.s_noise.';
            end
        end
    end
    
    
    %Choose method of error checking
    switch cfg.error_check
        case 'max'
            check_error=@max_error; %Choosing the maximum among system's dimensions
        otherwise
            check_error=@norm_error; %Euclidean norm among system's dimensions
    end
    
    %Initialization
    %Time step and error statistics
    cfg.dt=cell(cfg.Nt,cfg.Ntr);
    cfg.er=cell(cfg.Nt,cfg.Ntr);
    cfg.erPerD=cell(cfg.Nt,cfg.Ntr);
    cfg.erh=cell(cfg.Nt,cfg.Ntr);
    cfg.err=cell(cfg.Nt,cfg.Ntr);
    
    %Set the current time step as the initial one
    cdt=cfg.dt;
    
    %The current error
    cerr=1;
    
    %If the system is nonautonomous
    if strcmpi(cfg.systemType,'nonauto')
        %Initialization of the parameters structure array
        pout=struct(Nt,cfg.Ntr);
    end
    
    for jj=1:cfg.Ntr; %Start of trial loop
        
        %Initialization of the current trial
        
        %Store the initial condition
        x(1,:,jj)=x0;
        
        %Current x is the initial condition
        cx=x0;
        
        %Current time
        ct = t(1);
        
        %Meter of total iterations
        tot_steps=0;
        
        
        for ii=2:cfg.Nt; %Start of time loop
            
            if (cfg.Nlags>0)
                [ tlag,xlag,dxdtlag,plag ] = set_lags(ii,jj);
            else %...if there are no lags...
                xlag = [];
                dxdtlag = [];
                tlag = [];
                plag = [];
            end
            
            %Meter or steps taken between t(ii-1) and t(ii)
            step=1;
            
            while ct<t(ii) %Start of t(ii-1)->t(ii) time step loop
                
                
                %Initialize error statistics for this time step
                cfg.erPerD{ii,jj}(step,:) = zeros(1,cfg.D);
                cfg.er{ii,jj}(step) = 0;
                cfg.ers{ii,jj,step} = [];
                cfg.erh{ii,jj}(step) = 0;
                
                
                %Meter of iterations for the specific step
                it=1;
                
                ACCEPT=0;
                while ~ACCEPT %Start of iteration loop
                    
                    %Calculate first derivative
                    if strcmpi(cfg.systemType,'auto')
                        [dxdt4, dxdt5, ~]=Integrator(ct,cx,tlag,xlag,dxdtlag,cdt,System,p,plag);
                    else
                        [dxdt4, dxdt5, np]=Integrator(ct,cx,tlag,xlag,dxdtlag,cdt,System,p,plag); %get also next parameter setting (np)
                    end
                    %Calculate error
                    %4th order increment:
                    cdx = dxdt4*cdt;
                    x4 = cx+dx;
                    %5th order increment:
                    x5 = cx+dxdt5*cdt;
                    
                    %Error as the absolute difference
                    dx45= abs(x5-x4);
                    
                    %Check error
                    perr=cerr; %the current error (cerr) becomes the previous error (perr)
                    [ACCEPT, cerr, cfg.erPerD{ii,jj}(step,:), scale]= check_error(dx45,cfg.atol,cfg.rtol,cx);
                    
                    if ~ACCEPT %if the error is not acceptable
                        
                        %Store the excessive error and increase the meter of
                        %failed attempts (excessive  error hits) for this
                        %iteration
                        cfg.ers{ii,jj,step}(end+1) = cerr;
                        cfg.erh{ii,jj}(step) = cfg.erh{ii,jj}(step)+1;
                        
                        %Update time step and iteration indexes
                        it=it+1;
                        if (it<=cfg.max_iters)
                            cdt=update_dt(cdt,cerr,perr,cfg.dt_S,cfg.dt_a,cfg.dt_b,cfg.dt_min_ratio,cfg.dt_max_ratio,cfg.ndt,cfg.mdt);
                        else
                            %keyboard
                            disp(['t = ',num2str(ct)])
                            disp(['step = ',num2str(step)])
                            disp(['it = ',num2str(it)])
                            disp(['dt = ',num2str(cdt)])
                            disp(['Error per dimension = ',num2str(cfg.er{ii,jj}(step,:))])
                            disp(['Error = ',num2str(cerr)])
                            disp(['Error hits for this iteration= ',num2str(cfg.erh{ii,jj}(step))])
                            disp(['Errors for this iteration = ',num2str(cfg.ers{ii,jj}(step,:))])
                            disp(['atol = ',num2str(cfg.atol)])
                            disp(['rtol = ',num2str(cfg.rtol)])
                            disp(['scale=atol+rtol*abs(x) = ',num2str(scale)])
                            error('System stayed too many times on the same iteration without the error being reduced below error tolerance');
                        end
                    end
                    
                end %End of iteration loop
                
                %Store the current error of this successful step
                cfg.er{ii,jj}(step) = cerr;
                
                %Store the cdt
                cfg.dt{ii,jj}(step) = cdt;
                
                %Calculate next time point
                %We use the increment of 4th order in the absence of noise
                if NOISE
                    cdx = cdx + noiseFun(ct,dxdt4,cfg.s_noise,cdt);
                    cx = cx + cdx;
                else
                    cx = x4;
                end
                
                %Update parameters for nonautonomous parameters
                if strcmpi(cfg.systemType,'nonauto')
                    cp=np;
                end
                
                %Update time step
                cdt=update_dt(cdt,cerr,perr,cfg.dt_S,cfg.dt_a,cfg.dt_b,cfg.dt_min_ratio,cfg.dt_max_ratio,cfg.ndt,cfg.mdt);
                
                %Update iteration indexes
                step=step+1;
                tot_steps = tot_steps+1;
                
                %Update time
                pt=ct; %current time becomes previous time
                ct=ct+cdt;
                %Make sure we don't jump over the next time point t(ii)
                if sign(ct-t(ii))==sign(cdt)
                    ct=t(ii);
                    cdt = t(ii)-pt;
                end
                
                %Check if we exeeded the maximum number of total
                %interations allowed
                if (tot_steps>=cfg.max_tot_steps)
                    disp(['t = ',num2str(ct)])
                    disp(['tot_iters = ',num2str(step)])
                    disp(['dt = ',num2str(cdt)])
                    error('Maximum number of total steps allowed is exceeded');
                end
                
                
            end %End of time step loop
            
            %Store this time point
            x(ii,:,jj)=cx;
            dxdt(ii-1,:,jj)=dxdt4;
            if NOISE
                dx(ii-1,:,jj)=cdx;
            end
            
            %...and its parameters for nonautonomous systems
            if strcmpi(cfg.systemType,'nonauto')
                pout(ii-1,jj)=cp;
            end
            
            
        end %End of time loop
        
        if (cfg.Nlags>0)
            [ tlag,xlag,dxdtlag,plag ] = set_lags(ii+1,jj);
        else %...if there are no lags...
            xlag = [];
            dxdtlag = [];
            tlag = [];
            plag = [];
        end
            
        %Calculate first derivative for the last time point
        if strcmpi(cfg.systemType,'auto')
            [dxdt(ii,:,jj), ~, ~]=Integrator(ct,cx,tlag,xlag,dxdtlag,cdt,System,p,plag);
        else
            [dxdt(ii,:,jj), ~, pout(ii,jj)]=Integrator(ct,cx,tlag,xlag,dxdtlag,cdt,System,p,plag);
        end
        if NOISE
            dx(ii,:,jj) = dxdt(ii,:,jj)*cdt + noiseFun(ct,dxdt(ii,:,jj),cfg.s_noise,cdt);
        end
        
    end %End of trial loop
    
end

%Get rid of unused dimensions
x=squeeze(x);
dxdt=squeeze(dxdt);
if NOISE
    dx=squeeze(dx);
else
    dx=[];
end

%Store initial conditions as well
cfg.x0=x0;





%Plotting
if plotting>0
    D=cfg.D;
    Ntr=cfg.Ntr;
    tt=repmat(t,1,Ntr);
    iters=length(t);
    for ii=1:D;
        xx=x(:,ii,:);
        xx=xx(:);
        minX(ii)=min(xx);
        maxX(ii)=max(xx);
        
        dxx=dxdt(:,ii,:);
        dxx=dxx(:);
        mindX(ii)=min(dxx);
        maxdX(ii)=max(dxx);
    end
    clear xx dxx;
    
    
    figure((plotting-1)*10+1)
    for i=1:D;
        subplot(D,1,i);
        axis([-0.1 t(iters)+0.1 minX(i)-0.1 maxX(i)+0.1]);hold on;grid on;
        xx=squeeze(x(:,i,:));
        plot(tt,xx);hold on;
        ylabel(['x_',num2str(i)]);
        if (i==D)
            xlabel('t');
        end
        hold off;
    end
    
    figure((plotting-1)*10+2)
    for i=1:D;
        subplot(D,1,i);
        axis([-0.1 t(iters)+0.1 mindX(i)-0.1 maxdX(i)+0.1]);hold on;grid on;
        dxx=squeeze(dxdt(:,i,:));
        plot(tt,dxx);hold on;
        ylabel(['dx_',num2str(i),'/dt']);
        if (i==D)
            xlabel('t');
        end
        hold off;
    end
    
    if (D==1)
        figure((plotting-1)*10+3)
        axis([minX-0.1 maxX+0.1 mindX-0.1 maxdX+0.1]);hold on;
        xx=squeeze(x(:,i,:));
        dxx=squeeze(dxdt(:,i,:));
        plot(xx,dxx);hold on;
        xlabel('x');ylabel('dx/dt');hold off;
    elseif (D==2)
        figure((plotting-1)*10+3)
        axis([minX(1)-0.1 maxX(1)+0.1 minX(2)-0.1 maxX(2)+0.1]);hold on;
        xx1=squeeze(x(:,1,:));
        xx2=squeeze(x(:,2,:));
        plot(xx1,xx2);hold on;
        xlabel('x');ylabel('y');hold off;
    elseif (D==3)
        figure((plotting-1)*10+3)
        axis([minX(1)-0.1 maxX(1)+0.1 minX(2)-0.1 maxX(2)+0.1 minX(3)-0.1 maxX(3)+0.1]);hold on;
        xx1=squeeze(x(:,1,:));
        xx2=squeeze(x(:,2,:));
        xx3=squeeze(x(:,3,:));
        plot3(xx1,xx2,xx3);hold on;
        xlabel('x');ylabel('y');zlabel('z');hold off;
    elseif (D==4)
        figure((plotting-1)*10+3)
        xx1=squeeze(x(:,1,:));
        xx2=squeeze(x(:,2,:));
        xx3=squeeze(x(:,3,:));
        xx4=squeeze(x(:,4,:));
        subplot(1,2,1);
        axis([minX(1)-0.1 maxX(1)+0.1 minX(2)-0.1 maxX(2)+0.1 minX(3)-0.1 maxX(3)+0.1]);hold on;
        plot3(xx1,xx2,xx3);hold on;
        xlabel('x_1');ylabel('x_2');zlabel('x_3');hold off;
        subplot(1,2,2);
        axis([minX(3)-0.1 maxX(3)+0.1 minX(4)-0.1 maxX(4)+0.1 minX(1)-0.1 maxX(1)+0.1]);hold on;
        plot3(xx3,xx4,xx1);hold on;
        xlabel('x_3');ylabel('x_4');zlabel('x_1');hold off;
    end
    
    if NOISE
        
        
        clear dxx;
        
        figure((plotting-1)*10+4)
        for i=1:D;
            subplot(D,1,i);
            dxx=squeeze(dx(:,i,:));
            mindX(i)=min(dxx(:));
            maxdX(i)=max(dxx(:));
            axis([-0.1 t(iters)+0.1 mindX(i)-0.1 maxdX(i)+0.1]);hold on;grid on;
            plot(tt,dxx);hold on;
            ylabel(['dx_',num2str(i),'/dt']);
            if (i==D)
                xlabel('t');
            end
            hold off;
        end
    end
end


function cfg=createCFG
cfg.method='RK4';
cfg.systemType = 'auto';
cfg.s_noise=0;
cfg.timestep='fix';





function [dxdt,pout]=Euler(t,x,tlag,ylag,dydtlag,~,f,p,plag)%dt
[dxdt, pout]=f(t,x,tlag,ylag,dydtlag,p,plag);



function [dxdt,pout]=RK4(t,x,tlag,ylag,dydtlag,dt,f,p,plag)
[K1, pout] = f(t,x,          tlag,ylag,dydtlag,p,plag);
[K2, ~]    = f(t,x + dt*K1/2,tlag,ylag,dydtlag,p,plag);
[K3, ~]    = f(t,x + dt*K2/2,tlag,ylag,dydtlag,p,plag);
[K4, ~]    = f(t,x + dt*K3  ,tlag,ylag,dydtlag,p,plag);

dxdt = (K1 + 2*K2 + 2*K3 + K4)/6;



function [dxdt4, dxdt5, pout]=RK45(t,x,tlag,ylag,dydtlag,dt,f,p,plag)

%Constants of the Runge-Kutta Algorithm
a  =[1/4 3/32  1932/2197  439/216     -8/27];
b  =[    9/32 -7200/2197   -8          2];
c  =[          7296/2197 3680/513  -3544/2565];
d  =[                    -845/4104  1859/4104];
e  =                                 -11/40;
r4 =[25/216  1408/2565  2197/4104  -1/5];
r5 =[16/135 6656/12825 28561/56430 -9/50 2/55];


[k0, pout] = f(t,x,tlag,ylag,dydtlag,p,plag);

[k1, ~]    = f(t,x + a(1)*k0*dt,tlag,ylag,dydtlag,p,plag);

[k2, ~]    = f(t,x + (a(2)*k0 + b(1)*k1)*dt,tlag,ylag,dydtlag,p,plag);

[k3, ~]    = f(t,x + (a(3)*k0 + b(2)*k1 + c(1)*k2)*dt,tlag,ylag,dydtlag,p,plag);

[k4, ~]    = f(t,x + (a(4)*k0 + b(3)*k1 + c(2)*k2 + d(1)*k3)*dt,tlag,ylag,dydtlag,p,plag);

[k5, ~]    = f(t,x + (a(5)*k0 + b(4)*k1 + c(3)*k2 + d(2)*k3 + e*k4)*dt,tlag,ylag,dydtlag,p,plag);

%RK 4th order
dxdt4 = r4(1)*k0 + r4(2)*k2 + r4(3)*k3 + r4(4)*k4;

%RK 5th order
dxdt5 = r5(1)*k0 + r5(2)*k2 + r5(3)*k3 + r5(4)*k4 + r5(5)*k5;





function [ tlag,xlag,dxdtlag,plag ] = set_lagsFIX(ii)

global cfg t x dxdt x0 dxdt0 p pout;

%...initialize...
xlag = nan(cfg.D,cfg.Ntr,cfg.Nlags);
dxdtlag = nan(cfg.D,cfg.Ntr,cfg.Nlags);

%...get past values...
% xx = permute(x(:,:,1:(ii-1)),[3,1,2]);
% if ii>3
%     dxdx = permute(dxdt(:,:,1:(ii-2)),[3,1,2]);
% end
% tt = t(1:(ii-1));

%...calculate lagged times...
tlag = t(ii-1)-cfg.lags;
%...times before the initial time are set to the initial
%time...
ttau =tlag;
ttau(tlag<=t(1)) = t(1);

%...for each lag...
for iL = 1:cfg.Nlags;
    
    %...if it is the initial time...
    if (ttau(iL)==t(1))
        %...set the initial values as history...
        xlag(:,:,iL) = x0;
        if ~isempty(dxdt0)
            dxdtlag(:,:,iL) = dxdt0;
        end
        
    else
        
        %...if there are at least 2 past x values...
        if ii>2
            %...interpolate for each trial...
            %             for iTr = 1:cfg.Ntr;
            %                 xlag(:,iTr,iL) =  interp1(tt,squeeze(xx(:,:,iTr)),ttau(iL))';
            %             end
            xlag(:,:,iL) =  squeeze(interp1(t(1:(ii-1)),permute(x(:,:,1:(ii-1)),[3,1,2]),ttau(iL)));
        else
            %...set the initial value...
            xlag(:,:,iL) = x0;
        end
        
        %...if there are at least 2 past dxdt values...
        if ii>3
            %...interpolate for each trial...
            %             for iTr = 1:cfg.Ntr;
            %                 dxdtlag(:,iTr,iL) =  interp1(tt(1:end-1),squeeze(dxdx(:,:,iTr)),ttau(iL))';
            %             end
            dxdtlag(:,:,iL) =  squeeze(interp1(t(1:(ii-2)),permute(dxdt(:,:,1:(ii-2)),[3,1,2]),ttau(iL)));
        else
            if ~isempty(dxdt0)
                dxdtlag(:,:,iL) = dxdt0;
            end
            
        end
        
        
    end
    
end

%...parameters:
%...if the system is autonomous...
if strcmpi(cfg.systemType,'auto')
    plag = []; %...there is no special lagged p...
else
    %...if there at least 2 past pout values...
    if (size(pout,1)>3)&&(ii>3)
        %...for each lag...
        for iL = 1:cfg.Nlags;
            %...interpolate to the closest time index...
            this_i = round(interp1(t(1:(ii-2)),1:(ii-2),ttau(iL)));
            %...set the corresponding p...
            plag(iL,:) = pout(this_i,:);
        end
    else
        %...set initial p...
        plag = repmat(p,[cfg.Nlags,1]);
    end
end


function [ tlag,xlag,dxdtlag,plag ] = set_lags(ii,iTr)

global cfg t x dxdt x0 dxdt0 p pout;

%...initialize...
xlag = nan(cfg.D,cfg.Nlags);
dxdtlag = nan(cfg.D,cfg.Nlags);

%...get past values...
% xx = permute(x(:,iTr,1:(ii-1)),[3,1,2]);
% if ii>3
%     dxdx = permute(dxdt(:,iTr,1:(ii-2)),[3,1,2]);
% end
% tt = t(1:(ii-1));

%...calculate lagged times...
tlag = t(ii-1)-cfg.lags;
%...times before the initial time are set to the initial
%time...
ttau(tlag<=t(1)) = t(1);

%...for each lag...
for iL = 1:cfg.Nlags;
    
    %...if it is the initial time...
    if (ttau(iL)==t(1))
        %...set the initial values as history...
        xlag(:,iL) = x0(:,iTr);
        if ~isempty(dxdt0)
            dxdtlag(:,iL) = dxdt0(:,iTr);
        else
            dxdtlag=[];
        end
    else
        
        %...if there are at least 2 past x values...
        if ii>2
            %...interpolate for each trial...
%                 xlag(:,iTr,iL) =  interp1(tt,squeeze(xx(:,:,iTr)),ttau(iL))';
              xlag(:,iL) =  squeeze(interp1(t(1:(ii-1)),squeeze(x(1:(ii-1),:,iTr)),ttau(iL)));
        else
            %...set the initial value...
            xlag(:,iL) = x0(:,iTr);
        end
        
        
        %...if there are at least 2 past dxdt values...
        if ii>3
            %...interpolate for each trial...
            %                 dxdtlag(:,iTr,iL) =  interp1(tt(1:end-1),squeeze(dxdx(:,:,iTr)),ttau(iL))';
            dxdtlag(:,iL) =  squeeze(interp1(t(1:(ii-2)),squeeze(dxdt(1:(ii-2),:,iTr)),ttau(iL)));
        else
            if ~isempty(dxdt0)
                dxdtlag(:,iL) = dxdt0(:,iTr);
            end
        end
        
    end
    
end

%...parameters:
%...if the system is autonomous...
if strcmpi(cfg.systemType,'auto')
    plag = []; %...there is no special lagged p...
else
    %...if there at least 2 past pout values...
    if (size(pout,1)>3)&&(ii>3)
        %...for each lag...
        for iL = 1:cfg.Nlags;
            %...interpolate to the closest time index...
            this_i = round(interp1(t(1:(ii-2)),1:(ii-2),ttau(iL)));
            %...set the corresponding p...
            plag(iL,:) = pout(this_i,:);
        end
    else
        %...set initial p...
        plag = repmat(p,[cfg.Nlags,1]);
    end
end


        
function dxn = add_noise(~,dxdt,s_noise,dt)%t   
dxn=sqrt(dt)*s_noise.*randn( size(dxdt) );



function dxn = multi_noise(~,dxdt,s_noise,dt)%t
dxn=sqrt(dt)*s_noise.*dxdt.*randn( size(dxdt) );
    
    


function [ACCEPT maxError errorPerD scale] = max_error(er,atol,rtol,x)
scale = atol+rtol.*abs(x);
errorPerD = er./scale;
maxError = max(errorPerD);
ACCEPT = (maxError<=1);

function [ACCEPT normError errorPerD scale] = norm_error(er,atol,rtol,x)
scale = atol+rtol.*abs(x);
errorPerD = er./scale;
normError = sqrt( mean( (errorPerD).^2 ));
ACCEPT = (normError<=1);



function  cdt=update_dt(cdt,cerr,perr,S,a,b,dt_min_ratio,dt_max_ratio,ndt,mdt)
ratio = S*(cerr^a)*(perr^b);
ratio = min([ratio dt_min_ratio]);
ratio = max([ratio dt_max_ratio]);
cdt = cdt*ratio;
cdt = min([cdt,ndt]);
cdt = max([cdt,mdt]);








