function [p, x] = DPisoHistProb(y,Nbins,plotting,compare)

[N, M] = size(y);

if nargin<2
    Nbins = ceil(sqrt(N));
end

%Quantiles
p = [0:1/Nbins:1];
q = quantile(y,p);

%Bin centers
x = ( q(2:end,:) + q(1:(end-1),:) )/2; 

%The probability distribution as the derivative of the cumulative one...
p = 1./diff(q);
%...normalized to total probability one...
for iM = 1:M;
    p(:,iM) = p(:,iM)/trapz(x(:,iM),p(:,iM));  %http://stackoverflow.com/questions/5320677/how-to-normalize-a-histogram-in-matlab
end

if nargin>2
    if plotting
        plot(x,p)
        if compare
            [p2, x2] = hist(y,Nbins);
            for iM = 1:M;
                p2(:,iM) = p2(:,iM)/trapz(x2,p2(:,iM));  %http://stackoverflow.com/questions/5320677/how-to-normalize-a-histogram-in-matlab
            end
            hold on;
            plot(x2,p2,'--')
            hold off;
        end
    end
end
