function [Outs, hdr, mrk] = DPtfVAfun(EEGData, hdr, mrk,p)



% 
%Inputs:

%- EEGData: either loaded from a file or taken from MATLAB workspace during
%a VA MATLAB transform
%-hdr: header structure
%-mrk: marker structure
%-p: a structure of parameters of the function as described below:

%   -transform: the time-frequency transform to be used
%       either 'gabor' or 'gaborWavelet'
%   -param: parameter gamma of gabor transform or c of gabor wavelet
%   (Morlet)
%   -convORfft: flag selecting calculation via FFT or convolution in time domain 
%            either 'fft' or 'conv' (much slower option)
%   -flimSteps: frequemncy limits and number of steps, a vector of 
%          [lower frequency,  higher frequency,  number of frequency, steps]
%   -fscale: scale of frequency vector, 'LINEAR' or 'LOGARITHMIC'
%   -outputs: a cell of strings defining the output forms { 'complex',  'amplitude','power', 'phase'};
%%optionally:
%   -codepath: a string with the path to the code

%example script of setting the properties:
% procfun=@(EEGData, hdr, mrk,p)DPtfVAfun(EEGData, hdr, mrk,p);
% procparam.transform = 'gabor';
% procparam.param = 2*sqrt(pi);
% procparam.convORfft = 'fft';
% procparam.flimSteps = [      1                 81                   41];
% procparam.fscale = 'LINEAR';
% procparam.outputs = { 'complex'};
% procparam.codePath = {'\\Mpib10\InterBrain\EEGlab_VM\Denis\Software\DPtoolbox\Time-Frequency'}

%Outputs:
%Outs: a cell of the outputs in the same order as in p.outputs
%hdr: a cell of hdr structures, one for each output
%mrk: mrk structure, probably modified




%construct the frequency vector
%and set the corresponding properties in the hdr structure:
hdr.Layers = p.flimSteps(3);
hdr.LayerLowerLimit = p.flimSteps(1);
hdr.LayerUpperLimit = p.flimSteps(2);
hdr.LayerFunction = p.fscale;
clear flimSteps fscale
iF = 1:hdr.Layers;
if strcmp(hdr.LayerFunction,'LINEAR')
    %construct vector
    f = hdr.LayerLowerLimit + (hdr.LayerUpperLimit-hdr.LayerLowerLimit)/(hdr.Layers-1)*(iF-1);
elseif strcmp(hdr.LayerFunction,'LOGARITHMIC')
    %calculate log limits
    logLL = log(hdr.LayerLowerLimit);
    logUL = log(hdr.LayerUpperLimit);
    %construct vector
    f = exp( logLL+ (logUL-logLL)/(hdr.Layers-1)*(iF-1) );
end
p.f=f;

disp('...calculating transform...')
%tic
%Calculate transform with complex coefficients output:
if strcmpi(p.convORfft,'fft')
    Fs = 1000000/hdr.SamplingInterval;
    [Coef param] = DPtfviaFFT(EEGData,f,Fs,p.transform,p.param);
elseif strcmpi(convORfft,'conv')
    [Coef param] = DPtfviaCONV(EEGData,EEGTime,f,p.transform,p.param);
end
clear EEGData;

%define the parameters of the transform:
hdr.Prop1 = '';
hdr.Prop2 = '';
hdr.Prop3 = '';
hdr.Prop4 = '';
hdr.Prop5 = ['single,BrainVision.MorletFactor,',num2str(param)];
%toc

%Calculate and save outputs:
Nout0 = length(p.outputs);
GoodOut=[1:Nout0];
Outs = cell(1,Nout0);
hdrOrig = hdr;
hdr=cell(1,Nout0);
for iO = 1:Nout0;
    
    %Calculate output:
    if strcmpi(p.outputs{iO},'complex')
        disp(['...calculating ',p.outputs{iO},'...'])
        Outs{iO}.name = 'Coef';
        Outs{iO}.(Outs{iO}.name)=Coef;
        Outs{iO}.absORpow='abs';
        Outs{iO}.p=p;
        hdr{iO}=hdrOrig;
        hdr{iO}.DataType =  'TIMEFREQUENCYDOMAIN_COMPLEX';
        
    elseif strcmpi(p.outputs{iO},'amplitude')
        disp(['...calculating ',p.outputs{iO},'...'])
        %tic
        Outs{iO}.name = 'Ampl';
        Outs{iO}.(Outs{iO}.name) = abs(Coef);
        Outs{iO}.absORpow='abs';
        Outs{iO}.p=p;
        hdr{iO}=hdrOrig;
        hdr{iO}.DataType =  'TIMEFREQUENCYDOMAIN';
        %toc
        
        
    elseif strcmpi(p.outputs{iO},'power')
        disp(['...calculating ',p.outputs{iO},'...'])
        %tic
        Outs{iO}.name = 'Pow';
        Outs{iO}.(Outs{iO}.name) = abs(Coef).^2;
        Outs{iO}.absORpow='pow';
        Outs{iO}.p=p;
        hdr{iO}=hdrOrig;
        hdr{iO}.DataType =  'TIMEFREQUENCYDOMAIN';
        %toc
        
        
    elseif strcmpi(p.outputs{iO},'phase')
        disp(['...calculating ',p.outputs{iO},'...'])
        %tic
        Outs{iO}.name = 'Ph';
        Outs{iO}.(Outs{iO}.name) = angle(Coef);
        Outs{iO}.absORpow='abs';
        Outs{iO}.p=p;
        hdr{iO}=hdrOrig;
        hdr{iO}.DataType =  'TIMEFREQUENCYDOMAIN';
        %toc
else
        warning(['Ignoring output ',p.outputs{iO},'. It is not a valid output option, i.e. one of {''complex'',''ampl'',''power'',''phase''}']);
        GoodOut(1) = [];
    end
    
end
clear Coef;
%Check whether there is any valid output option at all:
Nout=length(GoodOut);
if Nout<1
    error('No valid output options')
end

disp('Done!')


