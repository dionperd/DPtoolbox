function [c, x, params] = DPcalcStableCDF_1D(y,method,Nbins,plotting) %, xi

%Inputs:
%-y: signals in a [points x signals] form
%-method: 
%non-parametric:
%         'empirical' 
%         'quantiles' for quantiles
%         the name of a kernel to be used with ksdensity canbe added to the
%         end of these two choices

%parametric:
%         'ASD' for parametric Alpha-Stable distributions toolbox
%            http://math.bu.edu/people/mveillet/html/alphastablepub.html#14
%         'LSD' for parametric Lévy Stable Distribution toolbox
%            http://www.ismm.ac.cn/ismmlink/LSD/LSD.htm
%'q' for using quantiles or
%nothing in order to calculate at all sample points.
%sampling, default: 'q'


%Optional: 
%-Nbins: number of bins
%-plotting: >0 for plotting

%Outputs
%-p: the probability distribution, 
%    normalized such that integral_x(p(x))dx=1 (trapz(x,p)=1)
%-x: the bin center points where the distribution is calculated
%
% %Optionally
% %-xi: the bin edge points where the distribution is calculated


[N, M] = size(y);

if nargin<3
    Nbins = ceil(sqrt(N));
end


if length(method)>4
      
    
    if  length(method)<=9 %if no kernel smoothing is required
        
        
        if strfind(method,'quantiles')
            
            %Quantiles and bin edges:
            c=[0:1/(Nbins-1):1].';

            %Bin edges:
            x = quantile(y,c);
                       
        else
            x = nan(N+1,M);
            c = x;
            for iM = 1:M;
                thisN = sum(~isnan(unique(y(:,iM)))) +1;
                [c(1:thisN,iM), x(1:thisN,iM)] = ecdf(y(:,iM));
            end
            c = c(2:end,:);
            x = x(2:end,:);
        end
        
        
    else
        
        kernel = method(10:end);
        
 
        if strfind(method(1:9),'quantiles')
            %Quantiles and bin edges:
            q=[0:1/(Nbins-1):1].';
            x = quantile(y,q);
            
            for iM = 1:M;
                %Cumulative distribution estimation:
                [c(:,iM), x(:,iM)] = ksdensity(y(:,iM),x(:,iM),'npoints',Nbins,'kernel',kernel,'function','cdf');
            end
        else
            for iM = 1:M;
                %Cumulative distribution estimation:
                yU = unique(y(:,iM));
                yU(isnan(yU))=[];
                thisN = length(yU);
                [c(1:thisN,iM), x(1:thisN,iM)] = ksdensity(y(:,iM),yU,'npoints',Nbins,'kernel',kernel,'function','cdf');
            end
        end
        
    end
    
    if nargout>2
        params = [];
    end
    
    
else
    
    if strcmpi(method(1:3),'ASD')
    
        fitfun = @stblfit;
        cdffun = @stblcdf;
        
    elseif strcmpi(method(1:3),'LSD')
        fitfun = @levystblfit;
        cdffun = @levystblcdf;
        
    else
        error('No valid method')
    end
    
    %Fit
    for iM = 1:M;
        params(:,iM) = fitfun(y(~isnan(y(:,iM)),iM));
    end
    
    %Calculate
    if length(method)<4
        
        x = nan(N,M);
        c = x;
        
        for iM = 1:M;
            
            yU = unique(y(:,iM));

            yU(isnan(yU))=[];
            thisN = length(yU);
            
            %Bin centers:
            x(1:thisN,iM) = yU;
            
            %Distribution calculation:
            c(1:thisN,iM) = cdffun(yU,params(1,iM),params(2,iM),params(3,iM),params(4,iM));
            
        end
        
        
    else
        
        %Quantiles:
        q=[0:1/(Nbins-1):1].';
        x = quantile(y,q);
        
        for iM = 1:M;
            %Distribution calculation:
            c(:,iM) = cdffun(x(:,iM),params(1,iM),params(2,iM),params(3,iM),params(4,iM));
        end
    end
    
    
end


if nargin>4
    if plotting
        plot(x,c)
    end
end
