%function DPVAfunFOLDER(InputPath,outputs,procfun,procparam,channels,trials,winRootPath,macRootPath,codePath)




%This function/script deals with the input-output processes between VA raw files
%that are all in the same folder, and MATLAB and calls a specific process function to be executed


%----------------------------for use as a script---------------------------

%For memory reasons, better to use as a script, so that EEGData is not
%repeatedly loaded in the memory

%Inputs:
if ispc
InputPath = '\\MPIB10\LIP-Git2\Denis\Karl\K_Export'; %win
end
if ismac
InputPath = '/Volumes/LIP-Git2/Denis/Karl/K_Export'; %mac
end
%Set the windows root folder:
winRootPath = '\\MPIB10';
macRootPath = '/Volumes';
codePath = {'\\MPIB10\InterBrain\EEGlab_VM\Karl\VAcoupl\DPtoolbox\io'};

channels='all';



%outputs and outputs folders
%this is a Nout x 3 cell of strings that defines:
%-1 if output data type is 1D (i.e. time domain), or 2 if it is 2D 
%    (e.g. Time-Frequency domain) (1st column)
%-the full paths of the folders where they should be saved (2nd column)
%-whether they should be saved in VA or MATLAB format (3rd column)
%2nd column: strings of folders full paths
%default=pwd
%3rd column:
%strings of either 'MATLAB' or 'VA' flags, default='MATLAB'
RawFolderPath=InputPath;
outputs = {...
    2,     RawFolderPath, 'VA';...
    };


%----------------------------for use as a script---------------------------

% % Add the DP io toolbox path
for ii=1:length(codePath);
    if ispc
        codePath{ii} = strrep(codePath{ii}, macRootPath, winRootPath);
%         InputPath = strrep(InputPath, macRootPath, winRootPath);
%         InputPath = strrep(InputPath, '/', '\');
        codePath{ii} = strrep(codePath{ii}, '/', '\');
    end
    if ismac
        codePath{ii} = strrep(codePath{ii}, winRootPath,macRootPath);
        codePath{ii} = strrep(codePath{ii}, '\', '/');
    end
    if isunix
        codePath{ii} = strrep(codePath{ii}, '\', '/');
    end
     addpath(codePath{ii})
end

trials='all';
procfun=@(EEGData, hdr, mrk,p)DPtfVAfun(EEGData, hdr, mrk,p);
p.transform = 'gabor';
p.param = sqrt(2*pi);
p.convORfft = 'fft';
p.flimSteps = [      1                 80                   20];
p.fscale = 'LINEAR';
p.outputs = { 'complex'};
p.codePath = codePath;

cd(InputPath);

FolderDir = dir;

for iFile = 3:length(dir);
    if any(strcmpi(FolderDir(iFile).name(end-2:end), {'eeg','dat','bin'}))
        disp(['torture file ',FolderDir(iFile).name])
        DPVAfun(outputs,procfun,p,channels,trials,FolderDir(iFile).name,winRootPath,macRootPath);
    end
end