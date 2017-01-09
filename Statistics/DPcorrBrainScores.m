%function [PLSres, PLScfg] = DPcorrBrainScores(PLSres,PLScfg,Bs,BsNames,plotting)
plotting=1;
load('PLSres.mat')
load('Bs.mat');

PLScfg.Nbs = length(BsNames);


for iBs= 1:PLScfg.Nbs; %For each behavioral or other subject measure...
    
    for iG = 1:PLScfg.Ng; %...for each group...
        
        x = Bs{iG}(:,iBs); %... get the measure (1 value per subject in a column)...
        
        for iM=1:PLScfg.Nm; %...for each measure...
            
            for iC = 1:PLScfg.Nc; %...for each condition...
                
                for iLV = 1:(PLScfg.Nlv); %...for each LV...

                    y = PLSres.(PLScfg.Measures{iM}).BrainScores{iG,iC}(:,iLV); %...get the subject's Brain Scores...
                    
                    %...calculate the correlation coefficient, the corresponding p-value...
                    %and the regression line:
                    [r, p, pol, R2] = DPlinRegress(x,y);
                    
                    %...and store them...
                    PLSres.(PLScfg.Measures{iM}).BrainScoresCorr{iG,iC,iLV}.(BsNames{iBs}).r = r;
                    PLSres.(PLScfg.Measures{iM}).BrainScoresCorr{iG,iC,iLV}.(BsNames{iBs}).p = p;
                    PLSres.(PLScfg.Measures{iM}).BrainScoresCorr{iG,iC,iLV}.(BsNames{iBs}).pol = pol;
                    PLSres.(PLScfg.Measures{iM}).BrainScoresCorr{iG,iC,iLV}.(BsNames{iBs}).R2 = R2;
                    clear y;
                end %iLV
            end %iC
            
        end %iG
    end%iM
    clear x;
end %iBs
PLScfg.Bs = Bs;
PLScfg.BsNames = BsNames;
save('PLSres.mat', 'PLSres', 'PLScfg');


if plotting
    
    Groups = {'Young','Old'};
    
    sbplot = [PLScfg.Ng*PLScfg.Nlv, PLScfg.Nc];
    
    for iM=1:PLScfg.Nm; %...for each measure...
        
        for iBs = 1:PLScfg.Nbs; %...for each behavioral or other subjects' measures...
            
            h=figure('Name',[PLScfg.Measures{iM},': Correlation and regression of Brain Scores with ',PLScfg.BsNames{iBs}]);
            
            for iLV = 1:(PLScfg.Nlv);
                
                for iG=1:PLScfg.Ng;
                    
                    x = PLScfg.Bs{iG}(:,iBs); %... get the measure (1 value per subject in a column)...
                    
                    for iC = 1:PLScfg.Nc;
                        
                        y = PLSres.(PLScfg.Measures{iM}).BrainScores{iG,iC}(:,iLV); %...get the subject's Brain Scores...
                        
                        r = PLSres.(PLScfg.Measures{iM}).BrainScoresCorr{iG,iC,iLV}.(PLScfg.BsNames{iBs}).r;
                        p = PLSres.(PLScfg.Measures{iM}).BrainScoresCorr{iG,iC,iLV}.(PLScfg.BsNames{iBs}).p;
                        pol = PLSres.(PLScfg.Measures{iM}).BrainScoresCorr{iG,iC,iLV}.(PLScfg.BsNames{iBs}).pol;
                        R2 = PLSres.(PLScfg.Measures{iM}).BrainScoresCorr{iG,iC,iLV}.(PLScfg.BsNames{iBs}).R2;
                        
                        yfit = polyval(pol,x);
                        
                        subplot(sbplot(1),sbplot(2),((iLV-1)*PLScfg.Ng + (iG-1))*sbplot(2)+iC)
                        plot(x,y,'o','Color',PLScfg.Colors{(iG-1)*PLScfg.Nc + iC});
                        hold on;
                        plot(x,yfit,'Color',PLScfg.Colors{(iG-1)*PLScfg.Nc + iC},'linewidth',1);
                        title(['LV ', num2str(iLV),...
                               ', ',Groups{iG},...
                               ', ', PLScfg.Conds{iC},...
                               ': r = ',num2str(r),...
                               ', p = ',num2str(p),...
                               ', R^2 = ',num2str(R2)]);
                        hold off;
                    end %iC
                end %iG
            end %iLV
            saveas(h,['BrScCorr_',PLScfg.Measures{iM},'_',PLScfg.BsNames{iBs},'.fig']);
            close(h);
            clear h;
        end %iBs
    end %iM
end %if
