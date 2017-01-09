function [p, x, params] = DPcalcStablePDF_1D(y,method,dq,Nbins,plotting) %, xi

%Inputs:
%-y: signals in a [points x signals] form
%-method: 
%non-parametric:
%         'equidist' for equidistant histogram
%         'isohisto' for iso-histogram
%         append the name of a kernel function, e.g. 'isohistonormal' for a
%         normal (gaussian) kernel, if kernel smoothing is required
%parametric:
%         'ASD' for parametric Alpha-Stable distributions toolbox
%            http://math.bu.edu/people/mveillet/html/alphastablepub.html#14
%         'LSD' for parametric Lévy Stable Distribution toolbox
%            http://www.ismm.ac.cn/ismmlink/LSD/LSD.htm
%add 'h' for sampling at equidistant points, 'q' for using quantiles or
%nothing in order to calculate at all sample points.
%sampling, default: 'q'
%-dq: quantile to be used in order to define the support of the
%    distribution, in the following way: q = [dq:(1-2*dq)/Nbins:1-dq]

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

if nargin<4
    Nbins = ceil(sqrt(N));
end


if length(method)>4
    
    
    if  length(method)<=8 %if no kernel smoothing is required
        
        
        if strfind(method,'equidist')
            
            q=[dq,1-dq];
            
            for iM = 1:M;
                %Boundaries:
                xi0 = quantile(y(:,iM),q);
                %                 %Bin edges:
                %                 xi(:,iM) = xi0(1):(xi0(2)-xi0(1))/Nbins:xi0(2);
                %Bin centers:
                %                 x(:,iM) = xi(1:(end-1),iM) + diff(xi(:,iM));
                temp = xi0(1):(xi0(2)-xi0(1))/(Nbins-1):xi0(2);
                if ~isempty(temp)
                    x(:,iM) = temp;
                    %Histogram:
                    [p(:,iM), x(:,iM)] = hist(y(:,iM),x(:,iM));
                    %Normalization:
                    p(:,iM) = p(:,iM)/trapz(x(:,iM),p(:,iM));  %http://stackoverflow.com/questions/5320677/how-to-normalize-a-histogram-in-matlab
                else
                    p(:,iM) = nan(Nbins,1);
                    x(:,iM) = zeros(Nbins,1);
                end
            end
            
            %
            %             if nargout>2
            %                 BinSiz = (x(2,iM)-x(1,iM));
            %                 BinSiz2 = BinSiz/2;
            %                 xi(:,iM) = [x(1,iM)-BinSiz2:BinSiz:x(end,iM)+BinSiz2].';
            %             end
            
            
            
        elseif strfind(method,'isohisto')
            
            %Quantiles and bin edges:
            q=[0,dq:(1-2*dq)/(Nbins-1):(1-dq),1].';
            
            %Bin edges:
            xi = quantile(y,q);
            
            
            %             % Euler method
            
            %             %Bin centers
            %             x = xi(1:(end-1),:) + diff(xi);
            %
            %             %The probability distribution as the derivative of the cumulative one...
            %             p = repmat(diff(q),1,M)./diff(xi);
            
            
            %Bin centers
            x = xi(2:(end-1),:);
            
            %The probability distribution as the derivative of the cumulative one...
            for iM=1:M;
                p(:,iM) = gradient(q(:,iM),xi(:,iM));
            end
            p = p(2:(end-1),:);
            p(isinf(p)) = nan;
            
        end
        
        
    else
        
        kernel = method(9:end);
        
        if strfind(method,'equidist')
            
            q=[dq,1-dq];
            
            for iM = 1:M;
                %Boundaries:
                xi0 = quantile(y(:,iM),q);
                %                 %Bin edges:
                %                 xi(:,iM) = xi0(1):(xi0(2)-xi0(1))/Nbins:xi0(2);
                %Bin centers:
                %                 x(:,iM) = xi(1:(end-1),iM) + diff(xi(:,iM));
                temp = xi0(1):(xi0(2)-xi0(1))/(Nbins-1):xi0(2);
                if ~isempty(temp)
                    x(:,iM) = temp;
                    %Density estimation:
                    [p(:,iM), x(:,iM)] = ksdensity(y(:,iM),x(:,iM),'kernel',kernel);
                else
                    p(:,iM) = nan(Nbins,1);
                    x(:,iM) = zeros(Nbins,1);
                end
            end
            %
            %             if nargout>2
            %                 BinSiz = (x(2,iM)-x(1,iM));
            %                 BinSiz2 = BinSiz/2;
            %                 xi(:,iM) = [x(1,iM)-BinSiz2:BinSiz:x(end,iM)+BinSiz2].';
            %             end
            
        elseif strfind(method,'isohisto')
            
            %Quantiles and bin edges:
            q=[0,dq:(1-2*dq)/(Nbins-1):(1-dq),1].';
            xi = quantile(y,q);
            
            for iM = 1:M;
                %Cumulative distribution estimation:
                [c(:,iM), xi(:,iM)] = ksdensity(y(:,iM),xi(:,iM),'npoints',Nbins,'kernel',kernel,'function','cdf');
            end
            
            
            %             % Euler method
            
            %             %Density calculation
            %             p = diff(c)./diff(xi);
            %             %Bin centers
            %             x = xi(1:(end-1),:) + diff(xi);
            
            %2nd order method
            %
            %Bin centers
            x = xi(2:(end-1),:);
            
            %The probability distribution as the derivative of the cumulative one...
            for iM=1:M;
                p(:,iM) = gradient(c(:,iM),xi(:,iM));
            end
            p = p(2:(end-1),:);
        end
        
    end
    
    if nargout>2
        params = [];
    end
    
    
else
    
    if strcmpi(method(1:3),'ASD')
        
        fitfun = @stblfit;
        pdffun = @stblpdf;
        
    elseif strcmpi(method(1:3),'LSD')
        fitfun = @levystblfit;
        pdffun = @levystblpdf;
        
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
        %         xi = x;
        p = x;
        
        for iM = 1:M;
            
            yU = unique(y(:,iM));
            
            yU(isnan(yU))=[];
            thisN = length(yU);
            %
            %             %Bin edges:
            %             xi(1:thisN,iM) = yU;
            %
            %Bin centers:
            x(1:thisN,iM) = yU;
            
            %Density calculation:
            if ~any(isnan(params(:,iM)))
                p(1:thisN,iM) = pdffun(yU,params(1,iM),params(2,iM),params(3,iM),params(4,iM));
            end
        end
        
    elseif strcmpi(method(4),'h')
        
        q=[dq,1-dq];
        
        for iM = 1:M;
            
            %             ymin = min(y(:,iM));
            %             ymax = max(y(:,iM));
            %             range = ymax - ymin;
            %             x(:,iM) = ymin:range/(Nbins-1):ymax;
            
            %Boundaries:
            xi0 = quantile(y(:,iM),q);
            
            %                 %Bin edges:
            %                 xi(:,iM) = xi0(1):(xi0(2)-xi0(1))/Nbins:xi0(2);
            %Bin centers:
            %                 x(:,iM) = xi(1:(end-1),iM) + diff(xi(:,iM));
            x(:,iM) = xi0(1):(xi0(2)-xi0(1))/(Nbins-1):xi0(2);
            
            %Density calculation:
            if ~any(isnan(params(:,iM)))
                p(:,iM) = pdffun(x(:,iM),params(1,iM),params(2,iM),params(3,iM),params(4,iM));
            end
            %             if nargout>2
            %                 BinSiz = (x(2,iM)-x(1,iM));
            %                 BinSiz2 = BinSiz/2;
            %                 xi(:,iM) = [x(1,iM)-BinSiz2:BinSiz:x(end,iM)+BinSiz2].';
            %             end
            
        end
        
    else
        
        %Quantiles:
        q=[dq:(1-2*dq)/(Nbins-1):(1-dq)].';
        x = quantile(y,q);
        p = nan(size(x));
        
        %Bin centers:
        %x(:,iM) = xi(1:(end-1),:) + diff(xi);
        %x = xi(2:(end-1),:);
        
        for iM = 1:M;
            %Density calculation:
            if ~any(isnan(params(:,iM)))
                p(:,iM) = pdffun(x(:,iM),params(1,iM),params(2,iM),params(3,iM),params(4,iM));
            end
        end
    end
    
    
end


if nargin>4
    if plotting
        plot(x,p)
    end
end




