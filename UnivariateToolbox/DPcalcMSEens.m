function [MSE, STD, G] = DPcalcMSEens(data,m,r,scales,num_tpts, num_trial, num_sc,downsample,normScale, slopeInds)  %override, DP excluded this input,
%        scales, DP excluded this output, 

% $$$ Written by Natasha Kovacevic - based on M. Costa's C code - simplified - skipped all sorts of command line options  - here we simultaneously get mse curves across trials or voxels
% $$$ Usage:   [mse, scales]= get_multiple_mse_curves_matlab(data,m,r,scales) 
% $$$ Inputs:
% $$$   data = time series data in (time variable) format, e.g (time trial) or (time voxel)
% $$$   m = pattern length (default 2)
% $$$   r = similarity criterion (default 0.5)
% $$$   scales = vector of scales (default [1:(floor(size(data,1)/50))])
% $$$   override = 0 or 1 (default is 0) - in some rare cases we want to oveeride this min50-rule
% $$$ Outputs:
% $$$   mse = entropy in (scale variable) format
% $$$   scales = return vector of scales actually used by the program

%-----------Modified by DP, find original in external folder---------------

%[num_tpts Ntr]= size(data);

%DP: to be moved in the parameters check of cfg prepare function

% if ~exist('m','var'), m=2; end
% if ~exist('r','var'), r=0.5; end
% if ~exist('scales','var'), scales = [1:(floor(num_tpts/50))] ; end

% % make sure that scales are positive integers
% scales = sort(unique(scales)); % sort scales in increasing order
% if sum(scales == round(scales)) ~= numel(scales)
%     error('scales vector must contain positive integers');
% end
% 
% if exist('override','var')
%     if ~override
%         % make sure that scales are  < (num_tpts/50)
%         scales = sort(intersect(scales,[1:(floor(size(data,1)/50))]));
%     end
% end


% normalize data
data=zscore(data);

MSE = zeros(num_sc,1); %(scale var)
STD = zeros(num_sc,1);
r0 = r;
for s=1:num_sc;
    sc = scales(s);
    
    % coarse grind time series at this scale
    num_cg_tpts = floor(num_tpts/sc);
    if downsample==1
        y = resample(data,num_cg_tpts,num_tpts,sc);
    else
        y = zeros(num_cg_tpts, num_trial);
        for t = 1:num_cg_tpts
            y(t,:) = mean(data((t-1)*sc + [1:sc],:),1);
        end
    end
    %DP: vectorize the result
    y=y(:); 
    
    STD(s) = std(y);
    if normScale
       r = r0*STD(s);
    end

    % calculate sample entropy of coarse ground time series y
    cont = zeros(m+1,1); %initialize counters
    
    %DP: slightly slower code that rejects patterns of transitions between trials 
    nlin_sc = num_cg_tpts - m; %number of m patterns per trial...
    nlin_scAlltr = num_trial*nlin_sc;%...and for all trials concatenated (excluding transition patterns)
    II0 = [1:nlin_sc].'; 
    II=zeros(nlin_scAlltr,1);
    II(1:nlin_sc)=II0;
    for iTr=2:num_trial;
        II((iTr-1)*nlin_sc+1 : iTr*nlin_sc) =(iTr-1)*nlin_sc+II0;
    end
        
    for ii = 1:nlin_scAlltr;
        i = II(ii);
        
        for ll=ii+1:nlin_scAlltr; % self-matches are not counted
            l=II(ll);
            
            k = 0;
            while ((k <= m) && (abs(y(i+k,1) - y(l+k,1)) <= r)) %DP changed <m to <=m...
                k = k + 1;
                cont(k,1) = cont(k,1) + 1;
            end
            %...and commented this useless part
            %             if ((k == m) && (abs(y(i+m,1) - y(l+m,1)) <= r))
            %                 cont(m+1,1) = cont(m+1,1) + 1;
            %             end
        end
    end

%     DP: Slightly faster code but it doesn't reject transition patterns between trials    
%     nlin_scAlltr = num_trial*num_cg_tpts-m;%...and for all trials concatenated (excluding transition patterns)
%     for i = 1:nlin_scAlltr;
%         
%         for l=i+1:nlin_scAlltr; % self-matches are not counted
%             
%             k = 0;
%             while ((k <= m) && (abs(y(i+k,1) - y(l+k,1)) <= r)) %DP changed <m to <=m...
%                 k = k + 1;
%                 cont(k,1) = cont(k,1) + 1;
%             end
%             %...and commented this useless part
%             %             if ((k == m) && (abs(y(i+m,1) - y(l+m,1)) <= r))
%             %                 cont(m+1,1) = cont(m+1,1) + 1;
%             %             end
%         end
%     end

    
    % calculate mse at this scale
    if (cont(m+1,1) == 0 || cont(m,1) == 0) %DP: why this???
        %DP: -log( 1 / ( (nlin_sc)*(nlin_sc -1) ) ) 
        %= log( nlin_sc*(nlin_sc-1) ) = log(nlin_sc^2- nlin_sc)
        MSE(s,1) = log(nlin_scAlltr^2 - nlin_scAlltr);
        %MSE(s,1) = inf;
    else
        MSE(s,1) = -log(cont(m+1,1)/cont(m,1));
    end
    
    
end 

if ~isempty(slopeInds)
    G = DPcalcSlope(log(scales),log(MSE),slopeInds);
else
    G=0;
end
