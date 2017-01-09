function [mse, mean_mse,std_mse,scales]= DGmse(data,m,r,scales) 
% $$$ Written by Natasha Kovacevic - based on M. Costa's C code - simplified - skipped all sorts of command line options
% $$$ Inputs:
% $$$   data = multiple trial data in (time trial) format
% $$$   m = pattern length (default 2)
% $$$   r = similarity criterion (defualt 0.5)
% $$$   scales = vector of scales (default [1:(floor(size(data,1)/50))])
% $$$ Outputs:
% $$$   mse: mse across scales for each trial (DP modification)
% $$$   mean_mse = mean mse across trials = vector of values, one value per scale
% $$$   std_mse = std mse across trials = vector of values, one value per scale
% $$$   scales = return vector of scales actually used by the program
%----------------------------Taken by Douglas Garrett----------------------
%---------------Modified by DP to calculate mse per trial------------------

  [num_tpts num_trials] = size(data);
  if ~exist('m','var'), m=2; end
  if ~exist('r','var'), r=0.5; end
  if ~exist('scales','var'), scales = [1:(floor(size(data,1)/50))] ; end
 
  % make sure that scales are positive integers < (num_tpts/50)
  scales = sort(unique(scales)); % sort scales in increasing order
  if sum(scales == round(scales)) ~= numel(scales)
    error('scales vector must contain positive integers');
  end
  scales = sort(intersect(scales,[1:(floor(size(data,1)/50))]));
  
  mse = zeros(numel(scales),num_trials);
	
  for trial=1:num_trials

    u = squeeze(data(:,trial)); % this trial time series 
    sd = std(u);
    r_new = r * sd;

    for s=1:numel(scales)

      sc = scales(s);
      % coarse grind time series at this scale
      num_cg_tpts = floor(num_tpts/sc);
      y = zeros(num_cg_tpts,1);
      for t=1:num_cg_tpts
        y(t) = mean(u((t-1)*sc + [1:sc]));
      end

      % calculate sample entropy of coarse ground time series
      nlin_sc = num_cg_tpts - m;
      cont = zeros(m+1,1);
      for i = 1:nlin_sc
        for l = (i+1):nlin_sc % self-matches are not counted
          k = 0;
          while k < m && abs(y(i+k) - y(l+k)) <= r_new
            k = k+1;
            cont(k) = cont(k) + 1;
          end
          if k == m && abs(y(i+m) - y(l+m)) <= r_new
            cont(m+1) = cont(m+1)+1;
          end
        end
      end
            
      % this is wasteful because we could get sample entropy for all m values [1:m],
			% instead I am only calculating for actual m
      if (cont(m+1) == 0 || cont(m) == 0)
        mse(s,trial) = -log(1/((nlin_sc)*(nlin_sc -1)));
      else
        mse(s,trial) = -log(cont(m+1)/cont(m));
      end
      

    end % for s=1:numel(scales)
  end % for trial=1:num_trials

  
  mean_mse = squeeze(mean(mse,2)); % mean entropy across trials
  std_mse = squeeze(std(mse,[],2)); % std entropy across trials
