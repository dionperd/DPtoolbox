function [MSE, RMSE, STD, MSEr, RMSEr,STDr, G, RG, Gr, RGr] = DPcalcMSEallEns(data,m,r,scales,num_tpts, num_trial, num_sc, slopeInds)  %override, DP excluded this input,
%        scales, DP excluded this output, 

% $$$ Written by Natasha Kovacevic - based on M. Costa's C code - simplified - skipped all sorts of command line options  - here we simultaneously get mse curves across trials or voxels
% $$$ Usage:   [mse, scales]= get_multiple_mse_curves_matlab(data,m,r,scales) 
% $$$ Inputs:
% $$$   data = time series data in (time 1iable) format, e.g (time trial) or (time voxel)
% $$$   m = pattern length (default 2)
% $$$   r = similarity criterion (default 0.5)
% $$$   scales = vector of scales (default [1:(floor(size(data,1)/50))])
% $$$   override = 0 or 1 (default is 0) - in some rare cases we want to oveeride this min50-rule
% $$$ Outputs:
% $$$   mse = entropy in (scale 1iable) format
% $$$   scales = return vector of scales actually used by the program

%-----------Modified by DP, find original in external folder---------------

%[num_tpts Ntr]= size(data);

%DP: to be moved in the parameters check of cfg prepare function

% if ~exist('m','1'), m=2; end
% if ~exist('r','1'), r=0.5; end
% if ~exist('scales','1'), scales = [1:(floor(num_tpts/50))] ; end

% % make sure that scales are positive integers
% scales = sort(unique(scales)); % sort scales in increasing order
% if sum(scales == round(scales)) ~= numel(scales)
%     error('scales vector must contain positive integers');
% end
% 
% if exist('override','1')
%     if ~override
%         % make sure that scales are  < (num_tpts/50)
%         scales = sort(intersect(scales,[1:(floor(size(data,1)/50))]));
%     end
% end


% normalize data
data=zscore(data);

MSE = zeros(num_sc,1); %(scale 1)
RMSE = zeros(num_sc,1); %(scale 1)
STD = zeros(num_sc,1);
MSEr = zeros(num_sc,1); %(scale 1)
RMSEr = zeros(num_sc,1); %(scale 1)
STDr = zeros(num_sc,1);

for s=1:num_sc;
    sc = scales(s);
    
    % coarse grind time series at this scale
    num_cg_tpts = floor(num_tpts/sc);
    y = zeros(num_cg_tpts, num_trial);
    for t = 1:num_cg_tpts
        y(t,:) = mean(data((t-1)*sc + [1:sc],:),1);
    end
    %resampling using a window length = sc
    yr = resample(data,num_cg_tpts,num_tpts,sc);
    
    %DP: vectorize the result
    y=y(:); 
    yr=yr(:);
    
    %Calculate STDs
    STD(s) = std(y);
    STDr(s) = std(y);

    %Adjust thresholds
    rSc = r*STD(s);
    rScR = r*STDr(s);
    
    % calculate sample entropy of coarse ground time series y
    cont = zeros(m+1,4); %initialize counters
    
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
            
            k = zeros(1,4);
           
            %MSE with coarse graining
            while ((k(1) <= m) && (abs(y(i+k(1),1) - y(l+k(1),1)) <= r0))
                k(1) = k(1) + 1;
                cont(k(1),1) = cont(k(1),1) + 1;
            end
            %MSE with resampling
            while ((k(2) <= m) && (abs(yr(i+k(2),1) - yr(l+k(2),1)) <= r0))
                k(2) = k(2) + 1;
                cont(k(2),2) = cont(k(2),2) + 1;
            end
            %RMSE with coarse graining
            while ((k(3) <= m) && (abs(y(i+k(3),1) - y(l+k(3),1)) <= rSc))
                k(3) = k(3) + 1;
                cont(k(3),3) = cont(k(3),3) + 1;
            end
            %RMSE with resampling
            while ((k(4) <= m) && (abs(yr(i+k(4),1) - yr(l+k(4),1)) <= rScR))
                k(4) = k(4) + 1;
                cont(k(4),4) = cont(k(4),4) + 1;
            end
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

    %MSE with coarse graining
    if (cont(m+1,1) == 0 || cont(m,1) == 0)
        %DP: -log( 1 / ( (nlin_sc)*(nlin_sc -1) ) )
        %= log( nlin_sc*(nlin_sc-1) ) = log(nlin_sc^2- nlin_sc)
        MSE(s,1) = log(nlin_scAlltr^2 - nlin_scAlltr);
    else
        MSE(s,1) = -log(cont(m+1,1)/cont(m,1));
    end
    %MSE with resampling
    if (cont(m+1,2) == 0 || cont(m,2) == 0)
        %DP: -log( 1 / ( (nlin_sc)*(nlin_sc -1) ) )
        %= log( nlin_sc*(nlin_sc-1) ) = log(nlin_sc^2- nlin_sc)
        MSEr(s,1) = log(nlin_scAlltr^2 - nlin_scAlltr);
    else
        MSEr(s,1) = -log(cont(m+1,2)/cont(m,2));
    end
    %RMSE with coarse graining
    if (cont(m+1,3) == 0 || cont(m,3) == 0)
        %DP: -log( 1 / ( (nlin_sc)*(nlin_sc -1) ) )
        %= log( nlin_sc*(nlin_sc-1) ) = log(nlin_sc^2- nlin_sc)
        RMSE(s,1) = log(nlin_scAlltr^2 - nlin_scAlltr);
    else
        RMSE(s,1) = -log(cont(m+1,3)/cont(m,3));
    end
    %RMSE with resampling
    if (cont(m+1,4) == 0 || cont(m,4) == 0)
        %DP: -log( 1 / ( (nlin_sc)*(nlin_sc -1) ) )
        %= log( nlin_sc*(nlin_sc-1) ) = log(nlin_sc^2- nlin_sc)
        RMSEr(s,1) = log(nlin_scAlltr^2 - nlin_scAlltr);
    else
        RMSEr(s,1) = -log(cont(m+1,4)/cont(m,4));
    end
    
end

if ~isempty(slopeInds)
        G = DPcalcSlope(log(scales),log(MSE),slopeInds);
        RG = DPcalcSlope(log(scales),log(RMSE),slopeInds);
        Gr = DPcalcSlope(log(scales),log(MSEr),slopeInds);
        RGr = DPcalcSlope(log(scales),log(RMSEr),slopeInds);
else
   G = 0;
   RG = 0;
   Gr = 0;
   RGr = 0;
end