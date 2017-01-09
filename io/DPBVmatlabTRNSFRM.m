%This function/script deals with the input-output processes between VA raw files
%and MATLAB when using VA MATLAB transform, 
%and calls a specific process function to be executed


%--------------------------for use as a function---------------------------
%function
%DPVAmatlabTRNSFRM(outputs,procfun,procparam,channels,trials,filename,NodeName,winRootPath,macRootPath,codePath,hdrFile,mrkFile) 
% %load workspace
% T = evalin('base','whos');
% for ii = 1:length(T)
%    C_ =  evalin('base',[T(ii).name ';']);
%    eval([T(ii).name,'=C_;']);
% end
% clear T C_ ii EEGData
%--------------------------for use as a function---------------------------



%----------------------------for use as a script---------------------------

%For memory reasons, better to use as a script, so that EEGData is not
%repeatedly loaded in the memory



%----------------------------for use as a script---------------------------

% %For memory reasons, better to use as a script, so that EEGData is not
% %repeatedly loaded in the memory
% 
% Inputs:

% %Set the windows root folder:
winRootPath = '\\Mpib10';
macRootPath = '/Volumes';
codePath = {'\\Mpib10\InterBrain\EEGlab_VM\Denis\Software\DPtoolbox\io'};

channels='all';
% 
trials='all';
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
outputs = {...
    2,     RawFolderPath, 'VA';...
    };
% 
%if there is a line in outputs as {1, '', 'return'}
%then that output is not saved but returned to VA

%hdrFile: 
%the full path FileName of the header file of the input dataset to be read. 
%This is optional since all the info can be read from the structure
%Properties of the VA MATLAB transform as it is done if hdrFile is not provided
%or if it is an empty string.
hdrFile='';

%mrkFile: 
%the full path filename of the marker file of the input dataset to be read. 
%This is optional since all the info can be read from the structure
%Markers of the VA MATLAB transform as it is done if mrkFile is not provided
%or if it is an empty string.

%if the user wants to provide a mrkFile but not a hdrFile, he has to
%provide an empty string for hdrFile
mrkFile='';


%set procfun and procparam

%----------------------------for use as a script---------------------------



%Add the DP io toolbox path
for ii=1:length(codePath);
    if ispc
        codePath{ii} = strrep(codePath{ii}, macRootPath, winRootPath);
    end
    if ispc
        codePath{ii} = strrep(codePath{ii}, '/', '\');
    end
    if ismac
        codePath{ii} = strrep(codePath{ii}, winRootPath,macRootPath);
    end
    if isunix
        codePath{ii} = strrep(codePath{ii}, '\', '/');
    end
    addpath(codePath{ii})
end

%Read the header and marker info:
disp('Reading data...')
%tic
[EEGDataNew,hdr, mrk, ~, ~, ~,~] = DPreadBVmatlabTRNSFRM(channels,trials,hdrFile,mrkFile);
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
[Outs hdr] = procfun(EEGDataNew, hdr, mrk,procparam);
%toc


%Save outputs:
Nout0 = size(Outs,1);

for iO = 1:Nout0;
    
    if strcmpi(outputs{iO,3},'return')
        
        if size(Outs{iO}.(Outs{iO}.name))==size(EEGData)
            EEGData = Outs{iO}.(Outs{iO}.name);
        end
            
    else
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
                DPwriteBVTF(single(Outs{iO}.(Outs{iO}.name)), fullfile(filepath,[FileName,'_',NodeName,'_',Outs{iO}.name,'.eeg']), hdr{iO}, Outs{iO}.absORpow, mrk);
            else
                %DPwriteBV(dat, filename, hdr, mrk)
                DPwriteBV(single(Outs{iO}.(Outs{iO}.name)), fullfile(filepath,[FileName,'_',NodeName,'_',Outs{iO}.name,'.eeg']), hdr{iO}, mrk);
            end
        else
            %save as MATLAB .mat file
            Res.(Outs{iO}.name) = Outs{iO}.(Outs{iO}.name);
            Res.p = Outs{iO}.p;
            Res.hdr = hdr{iO};
            Res.mrk = mrk;
            save(fullfile(filepath,[FileName,'_',NodeName,'_',Outs{iO}.name,'.mat']),'Res');
        end
    end
    %toc
end


disp('Done!')


