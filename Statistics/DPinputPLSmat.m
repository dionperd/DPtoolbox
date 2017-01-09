%function [PLSres, PLScfg] = DPrunPLStask( PLSfolder, RESfolder,  RESstruct, Conds, Groups, Measures, MeasOutsStruct, Ntr, Ns, Subjs)
%DPrunPLS This function calculates PLS analysis for several measures on a
%potentially multigroup and multisubject dataset

%   Outputs:
% -PLSres: structure with PLS results for each measure
% -PLScfg: structure with PLS configuration (PLSmat included)  for each measure

%% %---------------------A. Stack PLS matrix----------------------------------
%   Inputs:
RESfilePath = fullfile(pwd,'YA_OA_ori_vs_surr_analysis.mat');
%-RESfolder: full path to the file of the result structure
PLSmatFolder = pwd;%'full path to your PLS matrix folder';
StatsFolder = [PLSmatFolder,filesep,'Stats'];
if ~exist(StatsFolder,'dir')
    mkdir(StatsFolder)
end

load(RESfilePath)

RESstruct = 'C';
%-RESstruct: name of the result structure (string)
%-RESstruct: name of the result structure (string)
Conds={'Original','Surrogate'};%
%-Conds: cellstring of the names of conditions
Groups={'YA','OA'};
%-Groups: cellstring of the names of Groups 
%(a string common to all files, e.g., the extension, if there is only 1 group)
Ns= [31 28]; 
%-Ns: number of subjects of each group (vector of positive integers)
%load('C:\Users\CoordAgeEEG\Dropbox\Results\Subjects.mat')
Subjs =  {};%
SubjsFiles = {};
%-Subj: cell of cellstrings of the names of subjects for each group in case we want to use a
%       specific subset of them, otherwise Subjs={}
Measures = {'P','MRMSSD','MSE'};
%-Measures: a cellstring of the names of measures to be calculated
MeasOutsStruct = {[1 249], [1 100], [1 100]};
%-MeasOutsStruct: a cell of the same size as Measures with vectors
%                 corresponding to the size of the specific measure's 
%                 result PER TRIAL
Ntr=0;
%-Ntr: number of trials to be considered for the within subject average (integer),
%      if Ntr = 0, an average over all trials is taken

%---------------------For Plotting-----------------------------------------
%Inputs for plotting:
Colors = {rgb('DodgerBlue'),rgb('Red'), rgb('DarkBlue'),rgb('DarkRed')}; %rgb('ForestGreen'), rgb('Teal'),
LineStyle = '-';
varargin={}; %any other inputs you would like to give to the plotting function 
%plotMeasure = @(mean,std,Measure, Colors, Groups, Conds, Ns,varargin)plotMeanSTD_2D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
 plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns,LineStyle,varargin)plotMeanSTD_1D(mean,std,Measure, Colors, Groups, Conds, Ns,LineStyle,varargin);
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,Ns,varargin)plotMeanSTD_0D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
% % -plotMeasure: a handle to a function to plot each measure
%--------------------For Plotting------------------------------------------
%----------------------------End of Inputs---------------------------------


Nc = max(numel(Conds),1);
Ng = max(numel(Groups),1);
Nm =  numel(Measures);
Ncols = cellfun(@(x)prod(x),MeasOutsStruct);

Nr = Ns*Nc; %number of rows of PLSmat
%Initialize PLS matrices for each (sub)measure
for iM = 1:Nm;%For each measure...
    %Initialize
    PLSmat.(Measures{iM}) = cell(Ng,1);
    PLSmatStats.(Measures{iM}).mean = cell(Ng,Nc);
    PLSmatStats.(Measures{iM}).std = cell(Ng,Nc);
    for iG=1:Ng;  %...for each group...
        PLSmat.(Measures{iM}){iG} = reshape(C.(Measures{iM}){iG},Ns(iG)*Nc,Ncols(iM));
        for iC=1:Nc;
            %...calculate also statistics (mean and std) of all subjects for
            %this condition and group..
                RowIND = (iC-1)*Ns(iG)+[1:Ns(iG)];
                PLSmatStats.(Measures{iM}).mean{iG,iC} = mean( PLSmat.(Measures{iM}){iG}(RowIND,:) );
                PLSmatStats.(Measures{iM}).std{iG,iC} = std( PLSmat.(Measures{iM}){iG}(RowIND,:) );
                
                %...reshape to original sizes of the measures...
                PLSmatStats.(Measures{iM}).mean{iG,iC} = reshape( PLSmatStats.(Measures{iM}).mean{iG,iC}, MeasOutsStruct{iM} );
                PLSmatStats.(Measures{iM}).std{iG,iC} = reshape( PLSmatStats.(Measures{iM}).std{iG,iC}, MeasOutsStruct{iM} );
        end
    end
end


        
save([PLSmatFolder,filesep,'PLSmat.mat'],'PLSmat','PLSmatStats','Groups','Conds','Subjs','SubjsFiles','Measures','MeasOutsStruct',...
    'Ng','Nc','Ns','Nm','Ntr')
save([StatsFolder,filesep,'Stats.mat'],'PLSmatStats','Groups','Conds','Subjs','SubjsFiles','Measures','MeasOutsStruct',...
    'Ng','Nc','Ns','Nm','Ntr')
if exist('cfg','var')
    fields = fieldnames(cfg);
    for i=1:numel(fields);
        cfgMeasure.(fields{i}) = cfg.(fields{i});
    end
end
if exist('cfgMeasure','var')
    save([PLSmatFolder,filesep,'PLSmat.mat'],'cfgMeasure','-append')
end
if exist('Colors','var')
    save([PLSmatFolder,filesep,'PLSmat.mat'],'Colors','-append')
end

%---------------------For Plotting-----------------------------------------
for iM=1:Nm;
    h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM}, Colors, Groups, Conds, Ns,LineStyle,varargin);
    fileName = [StatsFolder, filesep, Measures{iM},'_stats.fig'];
    saveas(h,fileName);
    close(h);
    clear h;
end

for iM=1:Nm;
    h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM}, Colors, Groups, Conds, [],LineStyle,varargin);
    fileName = [StatsFolder, filesep, Measures{iM},'_stats.fig'];
    saveas(h,fileName);
    close(h);
    clear h;
end
%---------------------For Plotting-----------------------------------------

