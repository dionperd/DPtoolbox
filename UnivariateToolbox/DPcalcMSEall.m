function [MSE, RMSE, STD, MSEr, RMSEr,STDr, G, RG, Gr, RGr] = DPcalcMSEall(data,m,r,scales,num_tpts,num_vars,num_sc, slopeInds)  %override, DP excluded this input,
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
% 
% num_tpts = size(data,1);
% num_vars = size(data,2);

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



MSE = zeros(num_sc,num_vars); %(scale var)
MSEr = zeros(num_sc,num_vars); %(scale var)
RMSE = zeros(num_sc,num_vars); %(scale var)
RMSEr = zeros(num_sc,num_vars); %(scale var)
STD = zeros(num_sc,num_vars);
STDr = zeros(num_sc,num_vars);
for s=1:num_sc;
    sc = scales(s);
    
    % coarse grind time series at this scale
    num_cg_tpts = floor(num_tpts/sc);
    y = zeros(num_cg_tpts, num_vars);
    for t = 1:num_cg_tpts
        y(t,:) = mean(data((t-1)*sc + [1:sc],:),1);
    end
    %resampling using a window length = sc
    yr = resample(data,num_cg_tpts,num_tpts,sc);
    
    %calculate std for this scale
    STD(s,:) = std(y);
    STDr(s,:) = std(yr);
 
    r0 = r*STD(1,:); %threshold of original signal
    rSc = r*STD(s,:); %threshold adjusted to scale with coarse graining
    rScR = r*STDr(s,:); %threshold adjusted to scale with resampling
    
    
    % calculate sample entropy of coarse ground time series y
    nlin_sc = num_cg_tpts - m;
    cont = zeros(m+1,num_vars,4);
    for var = 1:num_vars
        for i = 1:nlin_sc
            for l = (i+1):nlin_sc % self-matches are not counted
                k = zeros(1,4);
                %MSE with coarse graining
                while ((k(1) <= m) && (abs(y(i+k(1),var) - y(l+k(1),var)) <= r0(var))) 
                    k(1) = k(1) + 1;
                    cont(k(1),var,1) = cont(k(1),var,1) + 1;
                end
                %MSE with resampling
                while ((k(2) <= m) && (abs(yr(i+k(2),var) - yr(l+k(2),var)) <= r0(var))) 
                    k(2) = k(2) + 1;
                    cont(k(2),var,2) = cont(k(2),var,2) + 1;
                end
                %RMSE with coarse graining
                while ((k(3) <= m) && (abs(y(i+k(3),var) - y(l+k(3),var)) <= rSc(var))) 
                    k(3) = k(3) + 1;
                    cont(k(3),var,3) = cont(k(3),var,3) + 1;
                end
                %RMSE with resampling
                while ((k(4) <= m) && (abs(yr(i+k(4),var) - yr(l+k(4),var)) <= rScR(var))) 
                    k(4) = k(4) + 1;
                    cont(k(4),var,4) = cont(k(4),var,4) + 1;
                end
            end
        end
    end
    
    % calculate mse at this scale
    for var = 1:num_vars;
        %MSE with coarse graining
        if (cont(m+1,var,1) == 0 || cont(m,var,1) == 0)
            %DP: -log( 1 / ( (nlin_sc)*(nlin_sc -1) ) )
            %= log( nlin_sc*(nlin_sc-1) ) = log(nlin_sc^2- nlin_sc)
            MSE(s,var) = log(nlin_sc^2 - nlin_sc);
        else
            MSE(s,var) = -log(cont(m+1,var,1)/cont(m,var,1));
        end
        %MSE with resampling
        if (cont(m+1,var,2) == 0 || cont(m,var,2) == 0)
            %DP: -log( 1 / ( (nlin_sc)*(nlin_sc -1) ) )
            %= log( nlin_sc*(nlin_sc-1) ) = log(nlin_sc^2- nlin_sc)
            MSEr(s,var) = log(nlin_sc^2 - nlin_sc);
        else
            MSEr(s,var) = -log(cont(m+1,var,2)/cont(m,var,2));
        end
        %RMSE with coarse graining
        if (cont(m+1,var,3) == 0 || cont(m,var,3) == 0)
            %DP: -log( 1 / ( (nlin_sc)*(nlin_sc -1) ) )
            %= log( nlin_sc*(nlin_sc-1) ) = log(nlin_sc^2- nlin_sc)
            RMSE(s,var) = log(nlin_sc^2 - nlin_sc);
        else
            RMSE(s,var) = -log(cont(m+1,var,3)/cont(m,var,3));
        end
        %RMSE with resampling
        if (cont(m+1,var,4) == 0 || cont(m,var,4) == 0)
            %DP: -log( 1 / ( (nlin_sc)*(nlin_sc -1) ) )
            %= log( nlin_sc*(nlin_sc-1) ) = log(nlin_sc^2- nlin_sc)
            RMSEr(s,var) = log(nlin_sc^2 - nlin_sc);
        else
            RMSEr(s,var) = -log(cont(m+1,var,4)/cont(m,var,4));
        end
    end
    
    
end % for s=1:numel(scales)

if ~isempty(slopeInds)
    for  var = 1:num_vars;
        G(var) = DPcalcSlope(log(scales),log(MSE(:,var)),slopeInds);
        RG(var) = DPcalcSlope(log(scales),log(RMSE(:,var)),slopeInds);
        Gr(var) = DPcalcSlope(log(scales),log(MSEr(:,var)),slopeInds);
        RGr(var) = DPcalcSlope(log(scales),log(RMSEr(:,var)),slopeInds);
    end
else
   G=zeros(1,num_vars);
   RG = G;
   Gr = G;
   RGr = G;
end