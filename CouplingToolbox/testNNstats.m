%import pastelsys.*;import pastelmath.*;import pastelgeometry.*;


pS.Nt=20;
pS.Ntr=2;
pS.D=2;
pS.X = (repmat( [1:pS.Nt].', [1,pS.D,pS.Ntr]) +0.001*randn(pS.Nt,pS.D,pS.Ntr) )/pS.Nt;
pS.X(:,1,2) = pS.X(:,1,2)+0.05; 
k=3;
r=0.002;
Wth=2;
maxDistanceSet=Inf;
treeMetric = 'euclidean';
queryIND = pS.Nt-3:pS.Nt;
L=length(queryIND);

figure(1)
trialColors = {'b','g','k','m','c'};
subplot(1,2,1)
for iTr = 1:pS.Ntr;    
    plot(pS.X(:,1,iTr),pS.X(:,2,iTr),[trialColors{iTr},'.']); hold on;
end
subplot(1,2,2)
for iTr = 1:pS.Ntr;    
    plot(pS.X(:,1,iTr),pS.X(:,2,iTr),[trialColors{iTr},'.']); hold on;
end

%TSTOOL
%Construct tree
tsS = DPconstrTREEtstool(pS,Wth,'euclidian');

%PASTEL
%Construct tree
plS = DPconstrTREEpastel(pS,Wth,treeMetric);


% %TSTOOL
% %Search N nearest neighbors
% [indNNtreeTS distNNts nnXts indNNtimeTS] = DPnnSearchTSTOOL(tsS, k, L, queryIND);
% 
% %PASTEL
% %Search N nearest neighbors
% [indNNtreePL distNNpl nnXpl indNNtimePL] = DPnnSearchPASTEL(plS, k, L, queryIND, treeMetric,maxDistanceSet);
% 
% 
% for iP = 1:L;
%     
%     iT = queryIND(iP);
%     
%     
%     for iTr = 1:pS.Ntr;    
%         
%         subplot(1,2,1)
%         h1(iP,iTr) = plot(pS.X(iT,1,iTr),pS.X(iT,2,iTr),[trialColors{iTr},'*'],'Markersize',10); hold on;
%         subplot(1,2,2)
%         h2(iP,iTr) = plot(pS.X(iT,1,iTr),pS.X(iT,2,iTr),[trialColors{iTr},'*'],'Markersize',10); hold on;
%         for iK = 1:k;
%             subplot(1,2,1)
%             hk1(iP,iTr,iK) = plot(nnXts(iP,1,iTr,iK),nnXts(iP,2,iTr,iK),'ro','Markersize',5); hold on;
%             axis([-0.1 1 -0.1 1.1])
%             subplot(1,2,2)
%             hk2(iP,iTr,iK) = plot(nnXpl(iP,1,iTr,iK),nnXpl(iP,2,iTr,iK),'ro','Markersize',5); hold on;
%             axis([-0.1 1 -0.1 1.1])
%             delete(hk1(iP,iTr,iK));
%             delete(hk2(iP,iTr,iK));
%         end
%         delete(h1(iP,iTr));
%         delete(h2(iP,iTr));
% 
%    end
% end


%TSTOOL
%Range search
[nncountTSrs indNNtreeTSrs distNNtsRS nnXtsRS indNNtimeTSrs] = DPrangeSearchTSTOOL(tsS, r, L, queryIND);

%PASTEL
%Range search
[nncountPL1] = DPrangeSearchPASTEL(plS, r, L, queryIND, treeMetric);
[nncountPL indNNtreePLrs distNNtsPLrs nnXplRS indNNtimePLrs] = DPrangeSearchPASTEL(plS, r, L, queryIND, treeMetric);



for iP = 1:L;
    
    iT = queryIND(iP);
    
    for iTr = 1:pS.Ntr;
        
        subplot(1,2,1)
        h1(iP,iTr) = plot(pS.X(iT,1,iTr),pS.X(iT,2,iTr),[trialColors{iTr},'*'],'Markersize',10); hold on;
        subplot(1,2,2)
        h2(iP,iTr) = plot(pS.X(iT,1,iTr),pS.X(iT,2,iTr),[trialColors{iTr},'*'],'Markersize',10); hold on;
        iT
        iTr
        if ( nncountTSrs(iP,iTr) == nncountPL(iP,iTr) )
            for iK = 1:nncountTSrs(iP,iTr);
                iK
                subplot(1,2,1)
                hk1{iP,iTr}(iK) = plot(nnXtsRS{iP,iTr}(iK,1),nnXtsRS{iP,iTr}(iK,2),'ro','Markersize',5); hold on;
                axis([-0.1 1.1 -0.1 1.1])
                
                subplot(1,2,2)
                hk2{iP,iTr}(iK) = plot(nnXplRS{iP,iTr}(iK,1),nnXplRS{iP,iTr}(iK,2),'ro','Markersize',5); hold on;
                axis([-0.1 1.1 -0.1 1.1])
                
                delete(hk1{iP,iTr}(iK));
                delete(hk2{iP,iTr}(iK));
            end
            nncountTSrs(iP,iTr)
            nncountPL(iP,iTr)
            nncountPL1(iP,iTr)
        else
            nncountTSrs(iP,iTr)
            nncountPL(iP,iTr)
            nncountPL1(iP,iTr)
            error('Shit!')
        end
        delete(h1(iP,iTr));
        delete(h2(iP,iTr));
        
    end
end