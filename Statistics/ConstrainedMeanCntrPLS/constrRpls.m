function [vbest, xbest, ibest, fbest, v, s,p,allPLSres,GA] = ...
    constrRpls(PLSfun, Nr,objfun, nonlcon,options)

global N;

Nrot = N*(N-1)/2;

lb = -90*ones(1,Nrot);
ub = -lb;
InitPop = zeros(1,Nrot);
options = gaoptimset('InitialPopulation',InitPop);
%x0 = zeros(1,Nrot);
for iR = 1:Nr;
    %[x(iR,:), fval(iR), exitflag(iR),output{iR}] = ...
        %patternsearch(objfun,x0,[],[],[],[],lb,ub, nonlcon,options);
    [GA.x(iR, :), GA.fval(iR), GA.exitflag(iR),GA.output{iR},...
        GA.population{iR}, GA.scores{iR}] = ...
        ga(objfun,Nrot,[],[],[],[],lb,ub, nonlcon,options);
    v(:,:,iR) = makeDesign(GA.x(iR,:));
    [s(:,iR), p(:,iR), allPLSres{iR}] = PLSfun(v(:,:,iR));
    save('C:\Users\dionperd\Dropbox\Dionysis\PostDocBerlin\StudentProjects\KarlGraphTheory\Denis\Optimization\GAprogress.mat','GA','v','allPLSres')
end

[fbest, ibest] = min(GA.fval);

xbest = GA.x(ibest,:)';
vbest = v(:,:,ibest);




