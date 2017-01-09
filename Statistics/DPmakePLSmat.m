%function [PLSres, PLScfg] = DPrunPLStask( PLSfolder, RESfolder,  RESstruct, Conds, Groups, Measures, MeasOutsStruct, Ntr, Ns, Subjs)
%DPrunPLS This function calculates PLS analysis for several measures on a
%potentially multigroup and multisubject dataset

%   Outputs:
% -PLSres: structure with PLS results for each measure
% -PLScfg: structure with PLS configuration (PLSmat included)  for each measure

%% %---------------------A. Stack PLS matrix----------------------------------
%   Inputs:
RESfolder='full path to your results folder';
%-RESfolder: full path to the folder of results .mat files
PLSmatFolder = 'full path to your PLS matrix folder';
StatsFolder = [PLSmatFolder,filesep,'Stats'];
if ~exist(StatsFolder,'dir')
    mkdir(StatsFolder)
end
%-PLSmatFolder: full path to the folder where PLS matrix will be saved
RESstruct = 'C';
%-RESstruct: name of the result structure (string)
Conds={'Cond1','Cond2','Cond3'};%
%-Conds: cellstring of the names of conditions
Groups={'Group1','Group2'};
%-Groups: cellstring of the names of Groups 
%(a string common to all files, e.g., the extension, if there is only 1 group)
Ns= [31 28]; 
%-Ns: number of subjects of each group (vector of positive integers)
%load('C:\Users\CoordAgeEEG\Dropbox\Results\Subjects.mat')
Subjs =  {};%
%-Subj: cell of cellstrings of the names of subjects for each group in case we want to use a
%       specific subset of them, otherwise Subjs={}
Measures = {'PLV','CCR'};
%-Measures: a cellstring of the names of measures to be calculated
MeasOutsStruct = {size(C.PLV), size(C.CCR)};
%-MeasOutsStruct: a cell of the same size as Measures with vectors
%                 corresponding to the size of the specific measure's 
%                 result PER TRIAL
Ntr=0;
%-Ntr: number of trials to be considered for the within subject average (integer),
%      if Ntr = 0, an average over all trials is taken

%---------------------For Plotting-----------------------------------------
%Inputs for plotting:
Colors = {rgb('DodgerBlue'),rgb('ForestGreen'),rgb('Red'), rgb('DarkBlue'),rgb('Teal'),rgb('DarkRed')};
varargin={}; %any other inputs you would like to give to the plotting function
plotMeasure = @(mean,std,Measure, Colors, Groups, Conds, Ns,varargin)plotMeanSTD_2D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,varargin)plotMeanSTD_1D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
% plotMeasure = @(mean,std,Measure, Colors, Groups, Conds,varargin)plotMeanSTD_0D(mean,std,Measure, Colors, Groups, Conds, Ns,varargin);
% % -plotMeasure: a handle to a function to plot each measure
%--------------------For Plotting------------------------------------------
%----------------------------End of Inputs---------------------------------




Nc = max(numel(Conds),1);
Ng = max(numel(Groups),1);
Nm =  numel(Measures);


Nr = Ns*Nc; %number of rows of PLSmat
%Initialize PLS matrices for each (sub)measure
for iM = 1:Nm;%For each measure...
    %Initialize
    PLSmat.(Measures{iM}) = cell(Ng,1);
    for iG=1:Ng;  %...for each group...
        PLSmat.(Measures{iM}){iG} = nan(Nr(iG),prod(MeasOutsStruct{iM}));
    end
    PLSmatStats.(Measures{iM}).mean = cell(Ng,Nc);
    PLSmatStats.(Measures{iM}).std = cell(Ng,Nc);
end

cd(RESfolder)

D = dir;
Nf=numel(D);
%The function to choose among the files...
if (isempty(Subjs))
    %...without subjects selection...
    ChooseFileCommand = '( ~isempty(strfind(D(iF).name,Groups{iG})) && ~isempty(strfind(D(iF).name,Conds{iC})) )';
else
    %...and with subjects selection...
    ChooseFileCommand = '( ~isempty(strfind(D(iF).name,Groups{iG})) && ~isempty(strfind(D(iF).name,Conds{iC}))&& ~isempty(strfind(D(iF).name,Subjs{iG}{iS})) )';
end

%A set of functions to get the data for each row of PLS mat:
if (Ntr==1) %this is the case where we already have 1 value for all trials (either mean or ensemble)
    for iM = 1:Nm;%...for each measure...
        %...get into temp the already averaged or ensemble value
        %for this subject...
        tempCommand{iM} = [RESstruct,'.',Measures{iM}];
    end
    
elseif (Ntr==0)%this is the case where we need to calculate the mean of ALL trials
    for iM = 1:Nm;%...for each measure...
        Mdim(iM) = length(MeasOutsStruct{iM});
        %...get into temp the average across all trials
        %for this subject...
        tempCommand{iM} = ['mean(',RESstruct,'.',Measures{iM},', Mdim(iM)+1)'];
    end
    
else
    for iM = 1:Nm;%...for each measure...
        Mdim(iM) = length(MeasOutsStruct{iM});
        %...construct the string that will index the
        %subset of selected trials...
        str{iM}='(';
        for iD = 1:Mdim(iM);
            str{iM} = [str{iM},':,'];
        end
        str=[str{iM},'[1:Ntr])'];
        %...get into temp the average across the first Ntr trials
        %for this subject...
        tempCommand{iM} = ['mean(',RESstruct,'.',Measures{iM},str{iM},', Mdim(iM)+1)'];
    end
end %if Ntr


tic
%Stack PLS matrix for each measure
disp('Stacking PLS matrix for each measure...')
for iG = 1:Ng; %For each group...
    disp('Group=')
    disp(Groups{iG})
    for iC = 1:Nc; %...for each condition...
        disp('Condition=')
        disp(Conds{iC})
        
        for iS=1:Ns(iG); %...for each subject...
            disp('Subject ')
            disp(iS)
            for iF=1:Nf; %...for each file...
                
                
                %...if this file corresponds to this group, condition (and subject)...
                if eval(ChooseFileCommand)

                    %...load the file...
                    load(D(iF).name)
                    
                    disp('PLSmat row index: ')
                    RowIND = (iC-1)*Ns(iG)+iS;
                    disp(RowIND)
                    
                    for iM = 1:Nm;%...for each measure...
                        
                        %...get this subject's data into temp...
                        temp=eval(tempCommand{iM});
                        
                        %                         disp(Measures{iM})
                        %                         disp(size(temp))
                        
                        %...finally add this subject's data into a row of PLS
                        %matrix:
                        ind = find(isnan(temp));
                        if ~isempty(ind)
                            warning('NaN values in PLS matrix converted to 0s')
                            iG
                            Groups{iG}
                            iC
                            Conds{iC}
                            iS
                            iF
                            D(iF).name
                            iM
                            Measures{iM}
                            temp(ind) = 0;
                        end
                        PLSmat.(Measures{iM}){iG}(RowIND,:) = temp(:);
                        
                    end  %for Nm
                    
                    %...remove the used file...
                    SubjsFiles{1,iG}{RowIND,1}=D(iF).name;
                    D(iF)=[];
                    Nf=numel(D);
                    break
                end% if choose file
            end %for Nf
        end%for Ns
        
        %...calculate also statistics (mean and std) of all subjects for
        %this condition and group..
        for iM = 1:Nm;%...for each measure...
            RowIND = (iC-1)*Ns(iG)+[1:Ns(iG)];
            PLSmatStats.(Measures{iM}).mean{iG,iC} = mean( PLSmat.(Measures{iM}){iG}(RowIND,:) );
            PLSmatStats.(Measures{iM}).std{iG,iC} = std( PLSmat.(Measures{iM}){iG}(RowIND,:) );
            
            %...reshape to original sizes of the measures...
            PLSmatStats.(Measures{iM}).mean{iG,iC} = reshape( PLSmatStats.(Measures{iM}).mean{iG,iC}, MeasOutsStruct{iM} );
            PLSmatStats.(Measures{iM}).std{iG,iC} = reshape( PLSmatStats.(Measures{iM}).std{iG,iC}, MeasOutsStruct{iM} );
            
        end
    end%for Nc
end %for Ng

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
    h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM}, Colors, Groups, Conds, Ns);
    fileName = [StatsFolder, filesep, Measures{iM},'_stats.fig'];
    saveas(h,fileName);
    close(h);
    clear h;
end

for iM=1:Nm;
    h = plotMeasure(PLSmatStats.(Measures{iM}).mean,PLSmatStats.(Measures{iM}).std,Measures{iM}, Colors, Groups, Conds, []);
    fileName = [StatsFolder, filesep, Measures{iM},'_stats.fig'];
    saveas(h,fileName);
    close(h);
    clear h;
end
%---------------------For Plotting-----------------------------------------

toc

clear all;
