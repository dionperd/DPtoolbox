function [s] = DPfcn_sl(data,params)
%FCN_SL        synchronization likelihood
%
%   S = FCN_SL(DATA,LAG,M,W1,W2,NREF);
%
%   Computes the synchronization likelihood (SL) for the recording series
%   in array, DATA (time x channel).
%
%   Inputs:     data,       recording series
%             params,       structure containing the following
%                lag,       lag parameter %DP: vector of absolute lags
%                  m,       embedding dimension
%                 w1,       starting window
%                 w2,       ending window
% %DP               nref,       number of recurrences in windows
%                 pref,     percentage of recurrences in windows
%
%   Outputs:       s,       SL matrix of dimensions time x channel pairings
%
%   Notes: this parameterization differs slightly from the traditional
%   parameterization of SL in that in does away with the (%DP: PREF) NREF parameter.
%   The parameters pref, nref, and w2 are all interrelated. Specifying any
%   pair of the two forces the remaining free parameter to be set.
%   Specifically, they are related in the following way:
%
%   (w2 - w1 - 1)*pref = nref;
%
%   This parameterization permits the function to perform somewhat faster.
%
%   References: Stam & van Dijk, (2002) doi: 10.1016/S0167-2789(01)00386-4.
%               Montez et al (2006) doi: 10.1016/j.neuroimage.2006.06.066.
%
%   Richard Betzel, Indiana University, 2012

%modification history
% 01.03.2012 - original 
%DP: 28-11-2014 - w2 equals all time series without boundary conditions
%               - embedding with variable lag per dimension

lag  = params.lag;
m    = params.m;
w1   = params.w1;
%DP--------------------
if mod(w1,2)~=0
    w1 = w1-1;
end
w12 = ceil(w1/2); 
%DP--------------------
w2   = params.w2;
if mod(w1,2)~=0
    w2 = w2-1;
end
w22 = ceil(w2/2); 

%nref = params.nref;
pref = params.pref;%DP

[npts,nchs]  = size(data);
%nptsadj      = npts - lag*(m - 1);
nptsadj      = npts - lag(end); %DP
%t            = 1:nptsadj;
%winsz        = (w2 - w1)*2 - 2;
%nsl          = length(w2:(nptsadj - w2 + 1));
masks        = find(triu(ones(nchs),1));
nmasks       = length(masks);

dummyadj         = zeros(nchs);
dummyx           = 1:(nchs^2);
dummyadj(dummyx) = dummyx;
dummyadj90       = rot90(dummyadj);
maskt            = dummyadj90(masks);

%s            = zeros(nsl,nmasks);
s            = zeros(nptsadj,nmasks);%DP
%embed        = zeros(m,nchs,nptsadj);
embed        = zeros(nptsadj,nchs,m);%DP

% for ipt = 1:nptsadj
%     embed(:,:,ipt) = data(ipt:lag:(ipt + lag*(m-1)),:);
% end
%DP--------------------
%...embed x per dimension
for iM = 1:m;
    embed(:,:,iM) = data( lag(iM) + 1 : lag(iM) + nptsadj,:);
end
%DP--------------------

clear data;

%embed = permute(embed,[3,2,1]);
%ind1  = ones(winsz,1);
ind2  = ones(1,1,nchs);
%nkeep = nref*2 + 1;

%count = 0;
%for i = w2:(nptsadj-w2+1);
 for i = 1:nptsadj;
    
    %count = count + 1;
    
%DP---------------------------------------------------------------------
    ind = [ [max(1,i-w22):(i-w12-1)] , [(i+w12+1):min(i+w22,nptsadj)]];%DP
    winsz  = length(ind);
    nref = pref*winsz;
    nrefo = floor(nref);
    if (nref~=nrefo)
        nref = nrefo;
        winsz0 = nref/pref;
        Dwinsz = winsz - winsz0;
        if (i==1)
            ind(end-Dwinsz+1:end)=[];
        else
            ind(1:Dwinsz)=[];
        end
        winsz = winsz0;
    end
    ind1  = ones(winsz,1);
    nkeep = nref + 1;
    
    eucd  = sum(bsxfun(@minus,embed(ind,:,:),embed(i,:,:)).^2,3).^0.5;
%DP---------------------------------------------------------------------
    %eucd  = sum(bsxfun(@minus,embed(abs(i-t) > w1 & abs(i-t) < w2,:,:),embed(i,:,:)).^2,3).^0.5;
    p     = sort(eucd);
    p     = p(nkeep,:);
    p     = p(ind1,:);
    h     = eucd - p < 0;
    h     = h(:,:,ind2);
    h1    = h(:,masks);
    h2    = h(:,maskt);
    
    %s(count,:) = sum(h1.*h2);
    s(i,:) = sum(h1.*h2)./nref; %DP
end
%s = s./(2*nref);