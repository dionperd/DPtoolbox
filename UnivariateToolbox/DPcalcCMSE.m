function [CMSE, CSTD, CG] = DPcalcCMSE(data,m,r,scales,num_tpts,num_vars,num_sc,normScale, slopeInds)  %override, DP excluded this input,
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


%-----------Modified by DP, find original in external folder---------------
%This function is modified to calculate composite MSE.

%Written by:
% Wu S-D, Wu C-W, Lin S-G, Wang C-C, Lee K-Y. 2013.
% Time Series Analysis Using Composite Multiscale Entropy.
% Entropy. 15:1069?1084.



CMSE = zeros(num_sc,num_vars); %(scale var)
CSTD = zeros(num_sc,num_vars); %(scale var)
for s=1:num_sc;
    
    sc = scales(s);
    
    for iSc = 1:sc; %for scale = sc we have sc possible ways to course grain
        
        % coarse grind time series at this scale
        %here we substract iSc-1 because we start course graining from point iSc
        num_cg_tpts = floor((num_tpts-iSc+1)/sc); 
        y = zeros(num_cg_tpts, num_vars);
        for t = 1:num_cg_tpts;
            %here we add iSc-1 because we start course graining from point
            %iSc
            y(t,:) = mean(data((t-1)*sc+iSc-1 + [1:sc],:),1); 
        end
        

        %calculate std for this scale and instance of course graining
        CSTDiSc(iSc,:) = std(y);
        
        if (s==1) %normalize threshold for the first scale's std (only one possible course graining)
            rSc = r*CSTDiSc(1,:);
        elseif normScale %and if we want the normalized version of MSE...
            %normalize for this scale (and course graining...) as well!
            rSc = r*CSTDiSc(iSc,:);
        end
        
        % calculate sample entropy of coarse ground time series y
        nlin_sc = num_cg_tpts - m;
        cont = zeros(m+1,num_vars);
        for var = 1:num_vars;
            for i = 1:nlin_sc;
                for l = (i+1):nlin_sc; % self-matches are not counted
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
                
                %here we devide by sc to take the average of all possible
                %sc ways to sample the data for scale sc
                CMSE(s,var) = CMSE(s,var)+log(nlin_sc^2 - nlin_sc)/sc; 
            else
                %here we devide by sc to take the average of all possible
                %sc ways to sample the data for scale sc
                CMSE(s,var) = CMSE(s,var) -log(cont(m+1,var)/cont(m,var))/sc;
            end
        end
        
    end
    
    CSTD(s,:)=mean(CSTDiSc);
    
end 


if ~isempty(slopeInds)
    for  var = 1:num_vars;
        CG(var) = DPcalcSlope(log(scales),log(CMSE(:,var)),slopeInds);
    end
else
   CG=zeros(1,num_vars);
end
