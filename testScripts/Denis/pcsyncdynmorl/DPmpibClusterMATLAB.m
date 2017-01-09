%Inputs:
username = 'perdikis';
processName = 'testCluster';
%Do:
ClusterPath = ['/home/mpib/',username]; %the path to the home folder in the cluster
clusterIP = [username,'@tardis'];
processClusterPath = [ClusterPath, '/',processName,'Folder'];%path to the process folder in your home folder in the cluster
system(['ssh ',clusterIP,' ''',ClusterPath,'; test -d ',processClusterPath,' || mkdir ',processClusterPath,'''']);
%system(['ssh ',clusterIP,' ''',ClusterPath,'; mkdir ',processClusterPath,'1''']);
resultsClusterPath = [processClusterPath, '/','Results'];%path to the results' folder in the cluster
system(['ssh ',clusterIP,' ''',ClusterPath,'; test -d ',resultsClusterPath,' || mkdir ',resultsClusterPath,'''']);

%IMPORTANT!: all paths should be given with the unix file separator: '/'
%In unix convention './' means current folder

%--------Copy .m file (and Toolboxes) and input data to the cluster--------
%Get the IP address of the PC...
%ifconfig | grep "inet " | grep -v 127.0.0.1 | cut -d\  -f2
[status result] = system('ifconfig | grep "inet " | grep -v 127.0.0.1 | cut -d\  -f2');
IP =  strtok(result);

%Create the .m file path and copy the file
%Inputs:
pcMfilePath ='/Users/dionysos/Dropbox/Dionysis/dptoolbox/testScripts/ClusterCode/testCluster.m'; %path to the .m file in your PC
%Do:
[pcMfileFolder, mFileName, mFileExt] = fileparts(pcMfilePath);
mFileClusterFolder = processClusterPath;%path to the .m file in your home folder in the cluster
mFileClusterPath = [mFileClusterFolder,'/',mFileName, mFileExt]; 
disp('.m file coping...')
%scp pcMfilePath username@tardis:mFileClusterPath
system(['scp ',pcMfilePath,' ',clusterIP,':',mFileClusterPath]);
disp('.m file copy DONE!')

%Create the data folder path and copy data folder
%Inputs:
pcDataPath ='/Users/dionysos/Dropbox/Dionysis/dptoolbox/testScripts/ClusterCode/pcDataPath'; %path to the data folder in your PC
%Do:
dataClusterPath = [processClusterPath,'/','inputData']; %the data path in the cluster
%Create input data path in the cluster if it is not already existing
system(['ssh ',clusterIP,' ''',ClusterPath,'; test -d ',dataClusterPath,' || mkdir ',dataClusterPath,'''']);
disp('Data folder coping...')
%scp -r  pcDataPath username@tardis:dataClusterPath
system(['scp -r ',pcDataPath,'/* ',clusterIP,':',dataClusterPath]);
disp('Data folder copy DONE!')

%Create the toolboxes' folder path and copy the folder
%Inputs:
pcToolboxesPath ='/Users/dionysos/Dropbox/Dionysis/dptoolbox/testScripts/ClusterCode/Toolboxes'; %path to the toolboxes' folder in your PC
if ~isempty(pcToolboxesPath)
    disp('Toolboxes'' folder coping...')
    toolboxesClusterFolder=[processClusterPath,'/Toolbox'];
    %Create toolbox path in the cluster if it is not already existing
    system(['ssh ',clusterIP,' ''',ClusterPath,'; test -d ',toolboxesClusterFolder,' || mkdir ',toolboxesClusterFolder,'''']);
    %scp -r  pcToolboxesPath username@tardis:pcToolboxesPath
    system(['scp -r ',pcToolboxesPath,'/* ',clusterIP,':',toolboxesClusterFolder]);
    disp('Toolboxes'' folder copy DONE!')    
    
end
%---------------Copy .m file input and data to the cluster-----------------



%-------------Create and send jobs to the cluster--------------------------
%Inputs:
Nj=3; %number of jobs
jobName = processName; %optional name of jobs
wallTime='24:00:00'; %hours, default value
dataSize = '1'; %gB

%Do:
%Create file path to the compiled sh script :
shFileName = ['jobs_',processName,'.sh'];
shPCFilePath = [pcMfileFolder,'/',shFileName];
shClusterFilePath = [processClusterPath,'/',shFileName];
shRunFileClusterPath = [ClusterPath,'/run_',mFileName,'.sh'];


%Open and write the script file
fid = fopen(shPCFilePath, 'w+');
%-----------------------------sh script------------------------------------ 
if ~isempty(pcToolboxesPath)
    %mcc -m mFileClusterPath -a toolboxesClusterPath
    fprintf(fid,['mcc -m %s -a %s\n'],mFileClusterPath,toolboxesClusterFolder);
else
    %mcc -m mFileClusterPath
    fprintf(fid,['mcc -m %s\n'],mFileClusterPath);
end

fprintf(fid,'for iJ in {1..%d}; do\n',Nj);
    
    %set the jobname variables:
    %echo "#PBS -N jobName$jobNumber" >> jobfile
    fprintf(fid,'echo "#PBS -N %s$iJ" >> jobfile\n',jobName);
    
    %set the resources of the jobs:
    %echo "#PBS -l walltime=30:00:00" >> jobfile
    fprintf(fid,'echo "#PBS -l walltime=%s" >> jobfile\n',wallTime);
    %echo "#PBS -l mem=50gb" >> jobfile
    fprintf(fid,'echo "#PBS -l mem=%sgb" >> jobfile\n',dataSize);
    
    %send the job to the cluster to run:
    %echo "#PBS -j oe" >> jobfile
    fprintf(fid,'echo "#PBS -j oe" >> jobfile\n');
    %echo "./run_KHVATF2CM_grid.sh /opt/matlab/interactive $ct" >> jobfile
    fprintf(fid,'echo "%s /opt/matlab/interactive $iJ" >> jobfile\n',shRunFileClusterPath);
    %qsub jobfile
    fprintf(fid,'qsub jobfile\n');
    %rm -f jobfile
    fprintf(fid,'rm -f jobfile\n');
        
fprintf(fid,['done']);
%-----------------------------sh script------------------------------------  
fclose(fid);
system(['scp ',shPCFilePath,' ',clusterIP,':',shClusterFilePath]);

%-------------------Connect to the cluster---------------------------------
%Do:
%ssh gridmaster2012 -l username
system(['ssh gridmaster2012 -l ',username]) %(substitute gridmaster2012 with 141.14.164.164 when DNS setup is broken)
%sh shRunFileClusterPath
%-------------------Connect to the cluster---------------------------------


%-------------Create and send jobs to the cluster--------------------------


%Exit
%exit %exit the cluster
%qstat %see what's going on with your jobs


%When the calculation is done, copy the output back to your PC 
%and free the space in the cluster! 
%Copy results:
pcResultsFolder =  '/Users/dionysos/Dropbox/Dionysis/dptoolbox/testScripts/ClusterCode/Results';
system(['test -d ',pcResultsFolder,' || mkdir ',pcResultsFolder])
system(['scp ',clusterIP,':',resultsClusterPath,'/* ',pcResultsFolder])
%Delte input data from cluster:
system(['ssh ',clusterIP,' ''',ClusterPath,'; rm -rf ',dataClusterPath,'''']);
%Delete results from cluster:
system(['ssh ',clusterIP,' ''',ClusterPath,'; rm -rf ',resultsClusterPath,'''']);


