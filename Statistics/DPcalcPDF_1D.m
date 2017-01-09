function [p, x] = DPcalcPDF_1D(y,method,Nbins,plotting) %, xi

%Inputs:
%-y: signals in a [points x signals] form
%-method: 'equidist' for equidistant histogram
%         'isohisto' for iso-histogram
%         append the name of a kernel function, e.g. 'isohistonormal' for a
%         normal (gaussian) kernel, if kernel smoothing is required
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

dq = 1/(Nbins-1);

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
            x(:,iM) = xi0(1):(xi0(2)-xi0(1))/(Nbins-1):xi0(2);
            %Histogram:
            [p(:,iM), x(:,iM)] = hist(y(:,iM),x(:,iM));
            %Normalization:
            p(:,iM) = p(:,iM)/trapz(x(:,iM),p(:,iM));  %http://stackoverflow.com/questions/5320677/how-to-normalize-a-histogram-in-matlab
        end
        %         if nargout>2
        %             BinSiz = (x(2,iM)-x(1,iM));
        %             BinSiz2 = BinSiz/2;
%             xi(:,iM) = [x(1,iM)-BinSiz2:BinSiz:x(end,iM)+BinSiz2].';
%         end
        
    elseif strfind(method,'isohisto')
        
        %Quantiles
        q=[0,dq:(1-2*dq)/(Nbins-1):(1-dq),1].';
        xi = quantile(y,q);
        
        
        %             % Euler method
        
        %             %Bin centers
        %             x = xi(1:(end-1),:) + diff(xi);
        %
        %             %The probability distribution as the derivative of the cumulative one...
        %             p = repmat(diff(q),1,M)./diff(xi);
        
        
        %2nd order method
        
        %Bin centers
        x = xi(2:(end-1),:);
        
        %The probability distribution as the derivative of the cumulative one...
        for iM=1:M;
            p(:,iM) = gradient(q(:,iM),xi(:,iM));
        end
        p = p(2:(end-1),:);

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
                x(:,iM) = xi0(1):(xi0(2)-xi0(1))/(Nbins-1):xi0(2);
                %Density estimation:
                [p(:,iM), x(:,iM)] = ksdensity(y(:,iM),x(:,iM),'kernel',kernel);
            end
        
%         if nargout>2
%             BinSiz = (x(2,iM)-x(1,iM));
%             BinSiz2 = BinSiz/2;
%             xi(:,iM) = [x(1,iM)-BinSiz2:BinSiz:x(end,iM)+BinSiz2].';
%         end
        
    elseif strfind(method,'isohisto')
        
        %Quantiles and bin edges:
        q=[0,dq:(1-2*dq)/(Nbins-1):(1-dq),1].';
        xi = quantile(y,q);
        
        for iM = 1:M;
            [c(:,iM), xi(:,iM)] = ksdensity(y(:,iM),xi(:,iM),'kernel',kernel,'function','cdf');
        end
        
        
        %             % Euler method
        
        %             %Density calculation
        %             p = diff(c)./diff(xi);
        %             %Bin centers
        %             x = xi(1:(end-1),:) + diff(xi);
        
        %2nd order method
        
        %Bin centers
        x = xi(2:(end-1),:);
        
        %The probability distribution as the derivative of the cumulative one...
        for iM=1:M;
            p(:,iM) = gradient(c(:,iM),xi(:,iM));
        end
        p = p(2:(end-1),:);
        %
    end
    
end



if nargin>3
    if plotting
        plot(x,p)
    end
end


function dxdy = diff2(x,y,M)

for iM=1:M;
    dxdy(:,iM) = gradient(y(:,iM),x(:,iM));
end

