function [PL nrm1 nrm2 dirin c1 c2 dirpar]=DPcalcPR(phi,fs,phaseTrnsfrm,method,ngrid,or)

% Rosenblum MG, Pikovsky A. 2001. 
% Detecting direction of coupling in interacting oscillators. 
% Physical Review E. 64:2?5.
% 
% Kralemann B, Cimponeriu L, Rosenblum MG, Pikovsky A, Mrowka R. 2008. 
% Phase dynamics of coupled oscillators reconstructed from data. 
% Physical Review E. 77:1?16.

if ~phaseTrnsfrm
    co_testproto(phi(:,1), phi(:,2));        % Testing protophases
    phi(:,1) = co_fbtransf1(phi(:,1));    % Preprocessing with univariate Fourier based transform
    phi(:,2) = co_fbtransf1(phi(:,2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%       WARNING,  if protophases are too strongly correlated
[dummy, PL n m]=co_maxsync(phi(:,1), phi(:,2), 2*or); 
if PL > 0.6;
    cprintf('Red','Synchronisation index is PLV>0.6: %f\n', num2str(PL));
    cprintf('Red','Protophases are strongly correlated, the results may not be reliable!\n');
    cprintf('Red','Calculation is stopped and nan values are returned!\n');
    nrm1 = nan;
    nrm2 = nan;
    dirin = nan;
    c1 = nan;
    c2 = nan; 
    dirpar= nan;
    return
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computing coeffs of Fourier expansion of the protophase dynamics
[Fcoef1,Fcoef2,f1,f2]=co_fexp2(phi(:,1),phi(:,2),or,fs,ngrid);


if strcmpi(method,'fourier')

    % DAMOCO Toolbox, function CO_FBTRANSF2, version 17.01.11
    % This high-level function performs all steps of the phase transformation at once,
    % using the Fourier-based method.
    % Given two protophases theta1 and theta2, the order of the Fourier expansion, or,
    % and the sampling frequency of the data, fs, this function computes:
    %      the true phases phi1 and phi2,
    %      the true coupling functions q1 and q2 on a grid of size ngrid,
    %      the Fourier coefficients of the coupling functions, Qcoef1 and Qcoef2,
    %      the frequencies (constant terms of coupling functions) omeg1 and omeg2,
    %      the norms of the coupling functions, norm1 and norm2,
    %      and the directionality index dirin.
    % Form of call:
    %      [phi1,phi2,omeg1,omeg2,q1,q2,Qcoef1,Qcoef2,norm1,norm2,dirin] = co_fbtransf2(theta1,theta2,or,fs,ngrid)
    %      [phi1,phi2,omeg1,omeg2,q1,q2,Qcoef1,Qcoef2,norm1,norm2,dirin] = co_fbtransf2(theta1,theta2,or,fs)
    % Input:
    %       theta1:     protophase of the 1st oscillator
    %       theta2:     protophase of the 2nd oscillator
    %       or:         order of the Fourier expansion: or = 10 is a good first choice.
    %                   To check and optimize this parameter: plot the  matrix of abs(Qcoef)
    %                   using co_plotcoef and verify that the coefficients at the boundary of the matrix are small.
    %                   If not, increase the parameter or. If the region, where the coefficients are nearly zero,
    %                   is large, reduce or.
    %       fs:       sampling frequency
    %       ngrid:      size of the grid to compute q1, q2; default value is ngrid=100
    % Output:
    %       phi1:       true phase of the 1st oscillator
    %       phi2:       true phase of the 2nd oscillator
    %       omeg1:      autonomous frequency (constant term of the coupling function) of the 1st oscillator
    %       omeg2:      autonomous frequency (constant term of the coupling function) of the 2nd oscillator
    %       q1:         coupling function of the 1st oscillator, computed on a grid
    %                   Note: It does not contain the constant term, i.e. the
    %                   phase dynamics is given by dphi1/dt = omeg1 + q1(phi1,phi2).
    %       q2:         coupling function of the 2nd oscillator, computed on a
    %                   grid, without constant term.
    %       Qcoef1:     the coefficients of the Fourier expansion of q1.
    %       Qcoef1:     the coefficients of the Fourier expansion of q2.
    %       norm1:      the norm of q1.
    %       norm2:      the norm of q2.
    %       dirin:      the index of directionality of interaction between both oscillators:
    %                   dirin = 1: unidirectional coupling from 1 to 2.
    %                   dirin = -1: unidirectional coupling from 2 to 1.
    %                   dirin = 0: equally strong bidirectional coupling
    %
    
    
   %%%%%%%%%%%%%%%%%%%%%           FOURIER BASED METHOD
   [sigfc1, sigfc2] = co_fbsolv(Fcoef1, Fcoef2);                       % Computing coefficients of phase transformation
   [phi(:,1), phi(:,2)] = co_fbth2phi(phi(:,1), phi(:,2), sigfc1, sigfc2);   % Performing phase transformation
   
   [Qcoef1, Qcoef2, q1, q2] = co_fexp2(phi(:,1), phi(:,2), or, fs, ngrid);   % Computing coupling functions and their Fourier coeffs
   omeg1 = real(Qcoef1(or+1,or+1));                    % Autonomous frequency is extracted from matrix Qcoef1
   omeg2 = real(Qcoef2(or+1,or+1));                    % Autonomous frequency is extracted from matrix Qcoef2
   q1 = q1 - omeg1;                                    % Coupling function of the phase does not contain the autonomous frequency
   q2 = q2 - omeg2;                                    % Coupling function of the phase does not contain the autonomous frequency
   Qcoef1(or+1,or+1) = 0;                              % The coefficients describe the coupling function only,
   Qcoef2(or+1,or+1) = 0;                              % they do not contain the autonomous frequency omeg1
   norm1 = co_fbnorm(Qcoef1);                          % Norms of the coupling functions
   norm2 = co_fbnorm(Qcoef2);
    

elseif strcmpi(method,'iter')
    
    % DAMOCO Toolbox, function CO_ITTRANSF2, version 17.01.11
    % This high-level function performs all steps of the phase transformation at once,
    % using the iteration method.
    % Given two protophases theta1 and theta2, the order of the Fourier expansion, or,
    % the sampling frequency of the data, fs,  and the grid size ngrid, this
    % function computes:
    %      the true phases phi1 and phi2,
    %      the true coupling functions q1 and q2 on a grid of size ngrid,
    %      the frequencies (constant terms of coupling functions) omeg1 and omeg2,
    %      the norms of the coupling functions, norm1 and norm2,
    %      and the directionality index dirin.
    % Form of call:
    %      [phi1,phi2,omeg1,omeg2,q1,q2,Nrmq1,Nrmq2,dirin] = co_ittransf2(theta1,theta2,or,fs,ngrid)
    %      [phi1,phi2,omeg1,omeg2,q1,q2,Nrmq1,Nrmq2,dirin] = co_ittransf2(theta1,theta2,or,fs)
    % Input:
    %       theta1:     protophase of the 1st oscillator
    %       theta2:     protopahse of the 2nd oscillator
    %       or:         order of the Fourier expansion: or = 10 is a good first choice.
    %       fs:      sampling frequency
    %       ngrid:      size of the grid to compute q1, q2; default value is ngrid=100
    % Output:
    %       phi1:       true phase of the 1st oscillator
    %       phi2:       true phase of the 2nd oscillator
    %       omeg1:      autonomous frequency (constant term of the coupling function) of the 1st oscillator
    %       omeg2:      autonomous frequency (constant term of the coupling function) of the 2nd oscillator
    %       q1:         coupling function of the 1st oscillator, computed on an a grid
    %                   Note: It does not contain the constant term, i.e. the
    %                   phase dynamics is given by dphi1/dt = omeg1 + q1(phi1,phi2).
    %       q2:         coupling function of the 2nd oscillator, computed on a
    %                   grid, without constant term.
    %       Nrmq1:      the norm of q1.
    %       Nrmq2:      the norm of q2.
    %       dirin:      the index of directionality of interaction between both oscillators:
    %                   dirin = 1: unidirectional coupling from 1 to 2.
    %                   dirin = -1: unidirectional coupling from 2 to 1.
    %                   dirin = 0: equally strong bidirectional coupling
    %
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%           Iteration Method
    niter=4;       % number of iterations; 4 is usually sufficient
    showiter=0;    % intermediate results are shown if showiter >0
    [q1,q2,omeg1,omeg2,norm1,norm2,sigma1,sigma2]=co_itersolv(f1,f2,niter,showiter);
    phi(:,1)=co_gth2phi(phi(:,1),sigma1);   % final transform of 1D-transformed protophases
    phi(:,2)=co_gth2phi(phi(:,2),sigma2);
    % q1=q1-omeg1; q2=q2-omeg2;            % coupling functions do not
    % contain constant terms

end

PL = co_sync(phi(:,1) ,phi(:,2),1,1);

% DAMOCO Toolbox, function CO_DIRIN, version 17.01.11
%
% Given the norms of the coupling functions of two coupled
% oscillators and their autonomous frequencies, this
% functions returns the directionality index.
% In case of symmetrical coupling dirin=0 holds, while purely
% uni-directional coupling yields dirin=1 or dirin=-1.
%
% Form of call: dirin = co_dirin(N1,N2,omeg1,omeg2)
% Input:        N1,N2 : norms of the coupling functions
%               omeg1,omeg2: frequencies
%
nrm1=norm1/omeg1; % strength of the external contribution to the phase
%DP: give nrm2 a negative value
nrm2=-norm2/omeg2; % dynamics normalized by the natural frequency
% to check for presence of interaction
S=nrm1-nrm2;
if S < 0.02
    cprintf('Red','nrm1+nrm2 < 0.02: %f\n', num2str(S));
    cprintf('Red','Warning: the coupling is very weak or the systems are not coupled!\n');
    cprintf('Red','Result on directionality index may be not reliable\n');
    cprintf('Red','Calculation is stopped and nan values are returned for PRn1, PRnn2, PRn!\n');
    nrm1 = nan;
    nrm2 = nan;
    dirin = nan;
else
    nrm1=nrm1/S;
    nrm2=nrm2/S;
    dirin= nrm1+nrm2; %  Directionality index
end

% DAMOCO Toolbox, function CO_DIRPAR, version 17.01.11
% Given the protophases theta1, theta2, this functions 
% returns the directionality index dirin, computed via 
% the partial derivatives of the coupling function with 
% respect to the external protophase.
% 
% Form of call:     dirin=co_dirpar(Fcoef1, Fcoef2)
%                  
% Input:
%       Fcoef1, Fcoef2:  Fourier coefficients for the model of phase
%                        dynamics for both systems
% Output:
%       dirin:           directionality index
%
S=size(Fcoef1); or=(S(1)-1)/2;
NP1=0;
NP2=0;
for n= -or : or;
    for m= -or : or;
        NP1 = NP1 + abs(1i*m*Fcoef1(n+or+1,m+or+1))^2;
        NP2 = NP2 + abs(1i*m*Fcoef2(n+or+1,m+or+1))^2;
    end;
end;
c1=sqrt(NP1);
c2=-sqrt(NP2); %DP: give c2 a negative value
S=c1 - c2;
if S < 0.02
    cprintf('Red','c1+c2 < 0.02: %f\n', num2str(S));
    cprintf('Red','Warning: the coupling is very weak or the systems are not coupled!\n');
    cprintf('Red','Result on directionality index may be not reliable\n');
    cprintf('Red','Calculation is stopped and nan values are returned for PRd1, PRd2, PRd!\n');
    c1 = nan;
    c2 = nan;
    dirpar = nan;
else
    c1 = c1/S;
    c2 = c2/S;
    dirpar= c1 + c2;
end


end