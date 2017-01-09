function [MSE, STD, G] = DPcalcMSE(data,m,r,scales,num_tpts,num_vars,num_sc,downsample, normScale, slopeInds)  %override, DP excluded this input,
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
STD = zeros(num_sc,num_vars);
for s=1:num_sc;
    sc = scales(s);
    
    % coarse grind time series at this scale
    num_cg_tpts = floor(num_tpts/sc);
    if downsample==1
        y = resample(data,num_cg_tpts,num_tpts,sc);
    else
        y = zeros(num_cg_tpts, num_vars);
        for t = 1:num_cg_tpts
            y(t,:) = mean(data((t-1)*sc + [1:sc],:),1);
        end
    end
    %calculate std for this scale
    STD(s,:) = std(y);
 
    if (s==1) %normalize threshold for the first scale's std
        rSc = r*STD(1,:);
    elseif normScale %and if we want the normalized version of MSE...
        %normalize for this scale as well!
        rSc = r*STD(s,:);
    end
    
    % calculate sample entropy of coarse ground time series y
    nlin_sc = num_cg_tpts - m;
    cont = zeros(m+1,num_vars);
    for var = 1:num_vars
        for i = 1:nlin_sc
            for l = (i+1):nlin_sc % self-matches are not counted
                k = 0;
                while ((k <= m) && (abs(y(i+k,var) - y(l+k,var)) <= rSc(var))) %DP changed <m to <=m...
                    k = k + 1;
                    cont(k,var) = cont(k,var) + 1;
                end
                %...and commented this useless part
%                 if ((k == m) && (abs(y(i+m,var) - y(l+m,var)) <= r))
%                     cont(m+1,var) = cont(m+1,var) + 1;
%                 end
            end
        end
    end
    
    % calculate mse at this scale
    for var = 1:num_vars;
        if (cont(m+1,var) == 0 || cont(m,var) == 0)
            %DP: -log( 1 / ( (nlin_sc)*(nlin_sc -1) ) )
            %= log( nlin_sc*(nlin_sc-1) ) = log(nlin_sc^2- nlin_sc)
            MSE(s,var) = log(nlin_sc^2 - nlin_sc);
        else
            MSE(s,var) = -log(cont(m+1,var)/cont(m,var));
        end
    end
    
end % for s=1:numel(scales)

if ~isempty(slopeInds)
    for  var = 1:num_vars;
        G(var) = DPcalcSlope(log(scales),log(MSE(:,var)),slopeInds);
    end
else
   G=zeros(1,num_vars);
end