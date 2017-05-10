function [f,causality,coherence,power] = compute_nparCausality(X,fs,fRes);
%Usage:  [f,causality,coherence,power] = compute_nparCausality(X,fs,fRes); 
%-------------------------------------------------------------------------------------------------
%Inputs: X = bivariate or multivariate signals in the form of 3D matrix  for time. trial. channel.  
%               fs = data sampling rate in Hz
%               fRes = a desired frequency resolution (e.g. 1)
%                            fRes is an optional input,  default value is fs/size(X,1)
%Outputs: f = frequencies at which causality (coherence or power) is computed
%                  causality = Granger causality between all pairs of signals (channels), e.g. 
%               :             causality(:,2,1) means causality from 2 to 1, and 1 to 2 is causality(:,1,2)
%               :             self-causality are set to zero, i.e. causality(:,k,k) = 0 for all k
%               : coherence = coherence between all pairs of signals (form: frequency. channel. channel)
%               : power      = 1-sided auto-spectra (form: frequency. channel)
%plot-example:   plot(f,causality(:,1,2),'b-',f,causality(:,2,1),'g--'); 
%-------------------------------------------------------------------------------------------------
% Ref : M. Dhamala, et al. , PRL (2007).
%-------------------------------------------------------------------------------------------------
%Flow of computation:  from signals (X) to 
%                     1. auto- & cross-spectra (S) 
%                     2. spectral matricial factors, then transfer function (H) & noise covariance (Z)
%                     3. Granger causality, coherence, 1-sided power (2*auto-spectra) 
%Written by M. Dhamala, USA, Aug 2006.
%-------------------------------------------------------------------------------------------------

[Nt, Ntr,Nc] = size(X); %Nt = number of timepoints, Ntr = trials, Nc = channels 

if nargin<3|fRes>fs/Nt,%a lower frequency resolution is achieved by zero-padding 
       fRes = fs/Nt; 
end 

[S,f]= sig2mTspect_nv(X,fs,fRes); %Not vectorized, less memory-demanding, signals to multitapered auto&cross-spectra 
%[S,f]= sig2mTspect(X,fs,fRes);% Vectorized routine, but memory-intensive, signals to multitapered auto&cross-spectra
%[S,f]= sig2spect(X,fs,fRes);%signals to auto-&cross-spectra without multitapering technique

if nargout>2
   spectra = permute(S,[3 1 2]);
   coherence = S2coh(spectra);
  if nargout>3
       for ichan = 1: Nc,
             power(:,ichan) = spectra(:,ichan,ichan);%one-sided power
       end
  end
end

if size(X,3)<3, % for pairwise causality of 2-channels
     [H,Z]=sfactorization_wilson(S,fs,f); %Wilson's algorithm for multivariate spectral matrix factorization
     % [H,Z]=sfactorization_bauer(S,fs,f); %Bauer's algorithm for multivariate spectral matrix factorization
     causality = hz2causality(H,S,Z,fs); 
else % for more than 2 channels, factorize two-channel spectra at a time and find pairwise causality 
     for ii = 1: Nc-1,
         for jj = ii+1: Nc, 
            S2 = S([ii jj],[ii jj],:);
            [H2,Z2]=sfactorization_wilson(S2,fs,f); %Wilson's algorithm for multivariate spectral matrix factorization
            % [H2,Z2]=sfactorization_bauer(S2,fs,f); %Bauer's algorithm for multivariate spectral matrix factorization
            cs = hz2causality(H2,S2,Z2,fs); 
            causality(:,ii,jj) = cs(:,1,2); causality(:,jj,ii) = cs(:,2,1);
         end
         causality(:,ii,ii) = 0; % self-causality is set to zero
     end 
     causality(:,Nc,Nc) = 0; %self-causality of the last channel is set to zero too   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------------------
% Functions below : sig2mTspect_nv.m, S2coh.m, sfactorization_wilson.m, hz2causality.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------------------

%----------------------------------sig2mTspect_nv.m------------------------------------------

function [S,f]= sig2mTspect_nv(X,fs,fRes);
%Usage: [S, f] = sig2mTspect_nv(X,fs,fRes);
%This function computes auto- & cross- spectra by using multitapers
%Inputs: X is multichannel data (a 3D-matrix in the form of time x trial x channel)
%               fs = sampling rate in Hz
%               nv stands for 'not vectorized program'
%               fRes = desired (lower) frequency resolution (e.g. 1), achieved with zero-padding 
%              default frequency resolution is fs/datalength
%Outputs: S = 3D matrix: m by m spectral matrix at each frequency point of f
%Note:     One can change nw (half the number of tapers) below and see the effect
%Written by M. Dhamala, USA, August 2006.

[N,Ntr,m] = size(X); % N = timepoints, Ntr = trials, m = channels 

if nargin<3|fRes>fs/N,
       npad = 0;  fRes = fs/N;  
end 
if (nargin==3) & (fRes<=fs/N), 
    npad = round((fs/fRes-N)/2);  %These many zeros will be padded on each side of the data
end
f = fs*(0:fix((N+2*npad)/2))/(N+2*npad);% upto Nyquist-f

nw = 7; % number of tapers = 2*nw-1 .......good nw are 1.5, 2, 3, 4, 5, 6, or 7...
[tapers,v] = dpss(N+2*npad,nw,2*nw-1);

S = zeros(m,m,N+2*npad);

for itrial = 1: Ntr,
     for ii = 1: m, Xft(:,:,ii) = mtfft(squeeze(X(:,itrial,ii)),tapers,fs,npad); end
     for ii = 1:m,
           for jj = 1:m,
                s(ii,jj,:) = squeeze(mean(Xft(:,:,ii).*conj(Xft(:,:,jj)),2));
                %averaging over tapers 
          end
     end
     S = S + s; 
end
S = S/Ntr; %averaging over trials
S = 2*S(:,:,1:fix(end/2)+1)/fs;%factor 2 for make one-sided spectra
S(:,:,1) = S(:,:,1)/2; %dc-power doesn't double for one-sided case
%--------------------------------------------------------------------------------------
function xf  = mtfft(data,tapers,fs,npad);
%Usage: xf = mtfft(data,tapers,fs,npad);
%Written by M. Dhamala (August 2006)

x0 = zeros(npad,size(data,2));
data = cat(1,x0,data); data= cat(1,data,x0);
data = data(:,ones(1,size(tapers,2)));
data = data.*tapers;xf = fft(data,[],1);

%---------------------------------------------------------------------------------

function [H, Z, psi] = sfactorization_wilson(S,fs,freq);
%Usage: [H, Z, psi] = sfactorization_wilson(S,fs,freq);
%Inputs: S (1-sided, 3D-spectral matrix in the form of Channel x Channel x frequency) 
%            : fs (sampling frequency in Hz)
%            : freq (a vector of frequencies) at which S is given
%Outputs: H (transfer function)
%       : Z (noise covariance)
%       : psi (left spectral factor)
%This function is an implemention of Wilson's algorithm (Eq. 3.1) for spectral matrix factorization
%Ref: G.T. Wilson,"The Factorization of Matricial Spectral Densities," SIAM J. Appl. Math.23,420-426(1972).
%Written by M. Dhamala & G. Rangrajan, USA, Aug 3-4, 2006.
%Email address: mdhamala@gsu.edu

m = size(S,1);N=length(freq)-1; tol = 1E-12; %tol is error-tolerence

%Step 1: Forming 2-sided spectral densities for ifft routine in matlab

f_ind=0;
for f=freq,
    f_ind=f_ind+1;
    Sarr(:,:,f_ind)=S(:,:,f_ind);
    if(f_ind>1),
        Sarr(:,:,2*N+2-f_ind)=S(:,:,f_ind).';
    end
end

% Step 2: Computing covariance matricies
for k1=1:m,
    for k2=1:m,
        gam(k1,k2,:)=real(ifft(squeeze(Sarr(k1,k2,:)))*fs);
    end
end

%Step 3: Initializing for iterations 
gam0 = gam(:,:,1);h = chol(gam0);

%%%%h = rand(m,m); h = triu(h); %arbitrary initial condition
 
for ind = 1: size(Sarr,3),
       psi(:,:,ind) = h; 
end

I = eye(m); % Defining m x m identity matrix

Niterations = 100; % Maximum number of iterations

% Step 4: Iterating to get spectral factors
 
for iter = 1: Niterations,
       
       for ind = 1: size(Sarr,3),
            g(:,:,ind)=inv(psi(:,:,ind))*Sarr(:,:,ind)*inv(psi(:,:,ind))'+I;%Eq 3.1
       end

       gp = PlusOperator(g,m,freq); %gp constitutes of positive and half of zero lags 

       psi_old=psi;
       for k = 1: size(Sarr,3),
             psi(:,:,k) = psi(:,:,k)*gp(:,:,k);
             psierr(k)=norm(psi(:,:,k)-psi_old(:,:,k),1);
       end
       psierrf=mean(psierr);   if(psierrf<tol),break;end; % checking convergence

end 

%for k = 1: length(freq),
%      Snew(:,:,k) = psi(:,:,k)*psi(:,:,k)'; % Snew: new spectral density
%end

%Step 5: Getting covariance matrix from spectral factors

for k1=1:m,
    for k2=1:m,
        gamtmp(k1,k2,:)=real(ifft(squeeze(psi(k1,k2,:))));
    end
end

% Step 6: Getting noise covariance & transfer function (see Example pp. 424)

A0=gamtmp(:,:,1); A0inv=inv(A0);

Z = A0*A0.'*fs; %Noise covariance matrix

for k = 1: length(freq),
      H(:,:,k) = psi(:,:,k)*A0inv; %Transfer function
end

%---------------------------------------------------------------------
function gp = PlusOperator(g,m,freq);
%This function is for [ ]+operation: 
%   to take the positive lags & half of the zero lag and reconstitute 
% M. Dhamala, August 2006

for k1=1:m,
    for k2=1:m,
          gam(k1,k2,:)= ifft(squeeze(g(k1,k2,:)));
    end
end

% taking only the positive lags and half of the zero lag

gamp = gam;beta0 = 0.5*gam(:,:,1); 
gamp(:,:,1) = triu(beta0);  %this is Stau
gamp(:,:,length(freq)+1:end) = 0;

% reconstituting
for k1=1:m,
    for k2=1:m,
         gp(k1,k2,:)= fft(squeeze(gamp(k1,k2,:)));
    end
end

%--------------------------------------------S2coh.m -----------------------

function coh = S2coh(S); 
%Input: S auto-& cross pectra in the form: frequency. channel. channel 
%Output: coh (Coherence) in the form: frequency. channel. channel 
%M. Dhamala, August 2006.

Nc = size(S,2);
for ii = 1: Nc,
   for jj = 1: Nc,
       coh(:,ii,jj) = real(abs(S(:,ii,jj)).^2./(S(:,ii,ii).*S(:,jj,jj)));
   end
end

%-------------------------------hz2causality.m ---------------------------- 

function causality = hz2causality(H,S,Z,fs);
%Usage: causality = hz2causality(H,S,Z,fs);
%Inputs: H = transfer function, S = 3-D spectral matrix;
%        Z = noise covariance,  fs = sampling rate
%Outputs: causality (Granger causality between all channels)
%               : auto-causality spectra are set to zero
% Reference: Brovelli, et. al., PNAS 101, 9849-9854 (2004).
%M. Dhamala, August 2006.

Nc = size(H,2);

for ii = 1: Nc,
    for jj = 1: Nc,
          if ii ~=jj,
              zc = Z(jj,jj) - Z(ii,jj)^2/Z(ii,ii);
              numer = abs(S(ii,ii,:));
              denom = abs(S(ii,ii,:)-zc*abs(H(ii,jj,:)).^2/fs);
              causality(jj,ii,:) = log(numer./denom);
          end
    end
    causality(ii,ii,:) = 0;%self-causality set to zero
end
causality = permute(causality,[3 1 2]); %freq x channel from x channel to

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% the subroutines below are alternative to some of the routines used %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------ sig2spect.m -------------------------------------------------------

function [S,f]= sig2spect(X,fs,fRes);
%Usage: [S, f] = sig2spect(X,fs,fRes);
%S = [Sxx Sxy; Syx Syy];Spectral matrix
%X is 3D-matrix: (time x trial x channel)
%Written by M. Dhamala, (June 2006).

[N,Ntr,m] = size(X); % N = timepoints, Ntr = trials, m = channels 

if nargin<3|fRes>fs/N,
       npad = 0;  fRes = fs/N;  
end 
if (nargin==3) & (fRes<=fs/N), 
    npad = round(fs/fRes-N);  %These many zeros will be padded on each side of the data
end
f = fs*(0:fix((N+npad)/2))/(N+npad);% upto Nyquist-f

Xft = fft(X,N+npad,1);

for ii = 1:m,
    for jj = 1:m,
        S(ii,jj,:) = mean(Xft(:,:,ii).*conj(Xft(:,:,jj)),2);%averaging over trials
    end
end
S = 2*S(:,:,1:fix(end/2)+1)/fs/N;%factor-2 to get one-sided
S(:,:,1) = 0.5*S(:,:,1);%factor-2 is not needed for dc

