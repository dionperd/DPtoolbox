function h = DPplot4parSearch(a,b,c,d,data,albl,blbl,clbl,dlbl,name,qlim,clickFun,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Na = length(a);
% Nb = length(b);
Nc = length(c);
Nd = length(d);

if nargin>11
   h = figure('Name',name,'units','normalized','outerposition',[0 0 1 1],'WindowButtonDownFcn', @plotThisPars);
else
   h = figure('Name',name,'units','normalized','outerposition',[0 0 1 1]);
end
iS = 0;
clim = quantile(data(:),qlim);
iM = floor(Nc/2);

for iC = 1:Nc;
    for iD = 1:Nd;
        iS = iS+1;
        subplot(Nc,Nd,iS)
        hold on;
        imagesc(a,b,squeeze(data(:,:,iC,iD)).',clim)
        axis tight
        %set(gca,'Ydir',normal)
        if all(iS~=[1:Nd:Nc*Nd])
            set(gca,'yticklabel',{})
            ylabel(blbl)
        else
            ylabel({blbl;[clbl,'=',num2str(c(iC))]})
        end
        if all(iS<=Nd*(Nc-1))
            set(gca,'xticklabel',{})
            xlabel(albl)
        else
            xlabel({[dlbl,'=',num2str(d(iD))];albl})
        end
        if iS==iM
            title(name)
        end
        if iS==Nc*Nd
            colorbar
        end
        set(gca,'userdata',[iC, iD])
        hold off
    end
end


    function plotThisPars(src, evt)
        
        ind = get(gca, 'CurrentPoint');
        thisA  = ind(1,1);
        thisB  = ind(1,2);
        ind = get(gca,'userdata');
        thisC = c(ind(1));
        thisD = c(ind(2));
        
        xlim = get(gca,'xlim');
        ylim = get(gca,'ylim');

        if ( (thisB>=ylim(1)) && (thisB<=ylim(2)) ) && ...
                ( (thisA>=xlim(1)) && (thisA<=xlim(2)) )
            
                [dum, ind] = min(abs(a-thisA));
                thisA = a(ind);
                [dum, ind] = min(abs(b-thisB));
                thisB = b(ind);
                clickFun(thisA,thisB,thisC,thisD,data,varargin);
        end
    end
end