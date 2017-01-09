function axP = plotBrainNet(net,axP,clim,axTitle,axXlbl,axYlbl,surfFile,nodeFile,optionsFile)

%It returns either the axis or the figure handle of the brain net

% nodeFile = '/media/dionysios/Windows/Users/dionperd/Dropbox/Dionysis/DPtoolbox/Plotting/external/BrainNetViewer/Data/ExampleFiles/Custom/EEG_21_interbrain.node';
% optionsFile ='/media/dionysios/Windows/Users/dionperd/Dropbox/Dionysis/DPtoolbox/Plotting/external/BrainNetViewer/Data/ExampleFiles/Custom/denisSymmetricInterbrainOptions.mat'; 
% surfFile = '';
%outputFile = '';

%save('net.edge','net');
%DPBrainNet_MapCfg('net.edge'nodeFile,,optionFile,surfFile); %H_BrainNet = 

DPBrainNet_MapCfg(net,nodeFile,optionsFile,surfFile); %H_BrainNet = 
hBN = gcf;

clear net;
if nargin>1
    hdummy = hBN;
    hBN = findall(hBN,'type','axes');
    %for iAx=1:2;
        copyobj(allchild(hBN(1)),axP); %hBN;
    %end
    clear hBN;
    close(hdummy);
    clear hdummy;
    axis(axP,'tight'); 
    set(axP,'clim',clim,'title',text('string',axTitle),'xlabel',text('string',axXlbl),'ylabel',text('string',axYlbl), 'visible', 'off');
    set(findall(axP, 'type', 'text'), 'visible', 'on');
    %colorbar;
    hold off; 
end

end