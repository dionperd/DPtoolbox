function DPVAfunFOLDER(fileIn)




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
OutputPath=[InputPath,'/MlOut']

%Set the windows root folder:
winRootPath = '\\MPIB10';
macRootPath = '/Volumes';
codePath = {'/Volumes/InterBrain/EEGlab_VM/Denis/Software/DPtoolbox/io'};

channels='all';%[57;59;58;41;49;43;45;47;22;24;26;28;30;6;8;10;12;14;1;2;3]
trials='all';

%%%%function to apply
procfun=@(EEGData, hdr, mrk,p)DPtfVAfun(EEGData, hdr, mrk,p);

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
        OutputPath = strrep(OutputPath, '\', '/');
        
    end
    if isunix
        codePath{ii} = strrep(codePath{ii}, '\', '/');
    end
    addpath(codePath{ii})
end

RawFolderPath=OutputPath;
outputs = {...
    2,     RawFolderPath, 'VA';...
    };

%%%prepare parameters
p.transform = 'gabor';
p.param = 2*sqrt(pi);
p.convORfft = 'fft';
p.flimSteps = [      2                 60                   30];
p.fscale = 'LINEAR';
p.outputs = { 'complex'};
p.codePath = {'\\MPIB10\InterBrain\EEGlab_VM\Denis\Software\DPtoolbox\Time-Frequency'};

cd(InputPath);

FolderDir = dir;

%%%%%use the subsequent code to evaluate single files only
ii=0;
for iFile = 3:length(dir);
    if any(strcmpi(FolderDir(iFile).name(end-2:end),{'eeg','dat','bin'}));%{'eeg','dat','bin'}
        ii=ii+1;
        files{ii}=FolderDir(iFile).name;
    end
end

for iFile=str2num(fileIn)
    disp('-----------------')
    disp(['TORTURE FILE ',files{iFile}])
    DPVAfun(outputs,procfun,p,channels,trials,files{iFile},winRootPath,macRootPath);
end


% %%%%use subsequent code to run all VA-file evaluations
% for iFile = 3:length(dir);
%     if any(strcmpi(FolderDir(iFile).name(end-2:end), {'dat'}));%{'eeg','dat','bin'}
%         disp('-----------------')
%         disp(['TORTURE FILE ',FolderDir(iFile).name])
%         DPVAfun(outputs,procfun,p,channels,trials,FolderDir(iFile).name,winRootPath,macRootPath);
%     end
% end


