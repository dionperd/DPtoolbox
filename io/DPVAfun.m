function DPVAfun(outputs,procfun,procparam,channels,trials,filename,winRootPath,macRootPath)

%This function/script deals with the input-output processes between VA raw files
%and MATLAB and calls a specific process function to be executed


%----------------------------for use as a script---------------------------

% %For memory reasons, better to use as a script, so that EEGData is not
% %repeatedly loaded in the memory
% 
% Inputs:

% filename = 'data.eeg';% or 'data.bin', or 'data.dat' 
% NodeName = 'Raw Data';
% %Set the windows root folder:
% winRootPath = '\\Mpib11';
% macRootPath = '/Volumes';
% codePath = {'\\Mpib10\InterBrain\EEGlab_VM\Denis\Software\DPtoolbox\io'};

% channels='all';
% 
% trials='all';
% 
% %outputs and outputs folders
% %this is a Nout x 3 cell of strings that defines:
% %-1 if output data type is 1D (i.e. time domain), or 2 if it is 2D 
% %    (e.g. Time-Frequency domain) (1st column)
% %-the full paths of the folders where they should be saved (2nd column)
% %-whether they should be saved in VA or MATLAB format (3rd column)
% %2nd column: strings of folders full paths
% %default=pwd
% %3rd column:
% %strings of either 'MATLAB' or 'VA' flags, default='MATLAB'
% outputs = {...
%     2,     RawFolderPath, 'VA';...
%     2,     RawFolderPath, 'VA';...
%     2,     RawFolderPath, 'VA';...
%     2,     RawFolderPath, 'VA';...
%     };
% 

%----------------------------for use as a script---------------------------


%Read the header and marker info:
disp('Reading data...')
%tic
[EEGData hdr mrk] = DPreadBV(filename, channels, trials);
%%%%ORIGINAL:::::[EEGData hdr mrk ~, ~, ~, ~] = DPreadBV(filename, channels, trials);
[~, FileName, ~] = fileparts(filename); 
%toc


%Calculate process
disp('...calculating process...')
%tic
%Add any more toolbox paths...
for ii=1:length(procparam.codePath);
    if isfield(procparam, 'codePath')
        if ispc
            procparam.codePath{ii} = strrep(procparam.codePath{ii}, macRootPath, winRootPath);
        end
        if ispc
            procparam.codePath{ii} = strrep(procparam.codePath{ii}, '/', '\');
        end
        if ismac
            procparam.codePath{ii} = strrep(procparam.codePath{ii}, winRootPath,macRootPath);
        end
        if isunix
            procparam.codePath{ii} = strrep(procparam.codePath{ii}, '\', '/');
        end
        addpath(procparam.codePath{ii})
    end
end
[Outs hdr] = procfun(EEGData, hdr, mrk,procparam);
%toc

%Save outputs:
Nout0 = length(Outs);

for iO = 1:Nout0;   
    
    %Define filepath...
    if exist(outputs{iO,2},'dir')
        filepath = outputs{iO,2};
        
    else
        filepath = pwd;
    end
    
    if ismac
        filepath = strrep(filepath, winRootPath, macRootPath);
    end
    if isunix
        filepath = strrep(filepath, '\', '/');
    end
    if ispc
        filepath = strrep(filepath,  macRootPath, winRootPath);
        filepath = strrep(filepath, '/','\');
    end
    %Save output
    disp(['...saving data of output...',num2str(iO)])
    %tic
    if strcmpi(outputs{iO,3},'VA')
        
        %write files in VA format:
        
        if outputs{iO,1}==2
            %DPwriteBVTF(dat, filename, hdr, absORpow, mrk)
            DPwriteBVTF(single(Outs{iO}.(Outs{iO}.name)), fullfile(filepath,[FileName,'_',Outs{iO}.name,'.eeg']), hdr{iO}, Outs{iO}.absORpow, mrk);
        else
            %DPwriteBV(dat, filename, hdr, mrk)
            DPwriteBV(single(Outs{iO}.(Outs{iO}.name)), fullfile(filepath,[FileName,'_',Outs{iO}.name,'.eeg']), hdr{iO}, mrk);
        end
    else
        %save as MATLAB .mat file
        Res.(Outs{iO}.name) = Outs{iO}.(Outs{iO}.name);
        Res.p = Outs{iO}.p;
        Res.hdr = hdr{iO};
        Res.mrk = mrk;
        
  
        save(([FileName,'_',Outs{iO}.name,'_.mat']),'Res');
    
        
    end
    %toc
end

disp('Done!')


