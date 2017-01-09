
%   Inputs:
%---------------------A. Stack PLS matrix----------------------------------
%RESfolder='C:\data\DATA_Aging_Complexity\ResultsFinal';
RESfolder='/media/dionysios/Windows/Users/dionperd/Dropbox/Dionysis/PostDocMRS/Results/Stats';
cd(RESfolder)
load('PLSmat.mat')
%-RESfolder: full path to the folder of results .mat files
Conds={'Rest EC','Odd notCount','Odd Count'};%'RestEO'
%-Conds: cellstring of the names of conditions
Groups={'Young','Old'};
%-Groups: cellstring of the names of Groups (empty if there are no groups)
%Measures = {'logP','PS','DOF','logV','VS','STD','MSE','logF','H'};
%Measures = {'logP','P','logV','STD','MSE','logF','H'};
Measures = {'H'};

%-Measures: a cellstring of the names of measures to be calculated
%           for instance: Measures = {'logP','PS','DOF','logV','VS','STD','MSE','logF','H'};
Nch=58; %number of channels here



Nc = numel(Conds);
Ng = numel(Groups);
Nm =  numel(Measures);


tic

%---------------------For Plotting-----------------------------------------
% Colors = {[69,139,116]/256,[72,61,139]/256,[255 64 64]/256,[255 165 79]/256, [127,255,212]/256,[100,149,237]/256,[165,42,42]/256,[255 127 36]/256};
%Colors = {[255 64 64]/256,[255 165 79]/256,[165,42,42]/256,[255 127 36]/256}; %TASK only
%Colors = {[69,139,116]/256,[72,61,139]/256,[127,255,212]/256,[100,149,237]/256}; %REST only
Colors = {rgb('DodgerBlue'),rgb('ForestGreen'),rgb('Red'), rgb('DarkBlue'),rgb('Teal'),rgb('DarkRed')};
%plotMeasure = @(x,Measure,cfgMeasure,Color,LineStyle,sbplt)plotUnivComplmultiPlotMeanSTD(x,Measure,cfgMeasure,Color,LineStyle,sbplt);
% plotMeasure = @(x,Measure,cfgMeasure,Color,LineStyle)plotUnivCompl(x,Measure,cfgMeasure,Color,LineStyle);
plotMeasure = @(mean,std,Measure,cfg,Color,LineStyle,sbplt)plotUnivComplmultiPlotMeanSTD(mean,std,Measure,cfg,Color,LineStyle,sbplt);
%-plotMeasure: a handle to a function to plot each measure
iX=0;
for iG=1:Ng;
    for iC=1:Nc;
        iX=iX+1;
        leg{iX} = [Groups{iG},'-',Conds{iC}];
    end
end
%---------------------For Plotting-----------------------------------------

for iM = 1:Nm;%...for each measure...
    
%     hm = figure('Name',[leg{iX},': ',Measures{iM},', subjects'' mean']);
%     hstd = figure('Name',[leg{iX},': ',Measures{iM},', subjects'' standard deviation']);
    h = figure('Name',[leg{iX},': ',Measures{iM},', subjects'' mean and standard deviation']);

    iX=0;
    %Plotting for each measure
    disp('Plotting for each measure...')
    for iG = 1:Ng; %For each group...
        disp('Group=')
        disp(Groups{iG})
        for iC = 1:Nc; %...for each condition...
            disp('Condition=')
            disp(Conds{iC})
   
            iX = iX+1;
                  
            %---------------------For Plotting-----------------------------------------
            %...and plot...
%             figure(hm)
%             hold on;
%             plotMeasure(PLSmatStats.(Measures{iM}).mean{iG,iC},Measures{iM},cfg,Colors{iX},'-',[2 4 iX]);
%             title(leg{iX})
%             hold off;
% 
%             figure(hstd)
%             hold on;
%             plotMeasure(PLSmatStats.(Measures{iM}).std{iG,iC},Measures{iM},cfg,Colors{iX},'-',[2 4 iX]);
%             title(leg{iX})
%             hold off;

              figure(h)
              hold on;
              plotMeasure(PLSmatStats.(Measures{iM}).mean{iG,iC},PLSmatStats.(Measures{iM}).std{iG,iC},Measures{iM},cfgMeasure,Colors{iX},'-',[1 Nc iC]); %[Ng Nc iX]
              %title(leg{iX})
              title([Conds{iC},', Cz'])
              legend(Groups)
              
              %hold off;
              
            %---------------------For Plotting-----------------------------------------
            
            
        end %for Nc
        
    end%for Ng
    
%     fileNameM = [Measures{iM},'_mean'];
%     saveas(hm,[fileNameM,'.fig']);
%     close(hm);
%
%     fileNameSTD = [Measures{iM},'_std'];
%     saveas(hstd,[fileNameSTD,'.fig']);
%     close(hstd);

    %fileName = [Measures{iM}];
    saveas(h,[Measures{iM},'_perCond.fig']);
    close(h);
end %for Nm
toc
