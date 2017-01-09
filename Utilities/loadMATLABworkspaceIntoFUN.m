%load workspace
T = evalin('base','whos');
for ii = 1:length(T)
   C_ =  evalin('base',[T(ii).name ';']);
   eval([T(ii).name,'=C_;']);
end
clear T C_ ii EEGData