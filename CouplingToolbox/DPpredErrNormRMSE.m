function Delta = DPpredErrNormRMSE(x,xPred)

%This function calculate the performance of a predictor via the normalized
%root mean squared prediction error

% Schiff SJ, So P, Chang T, Burke RE, Sauer T. 1996. 
% Detecting dynamical interdependence and generalized synchrony through 
% mutual prediction in a neural ensemble. 
% Phys Rev E Stat Phys Plasmas Fluids Relat Interdiscip Topics.
% 54:6708?6724.
% 
% Le Van Quyen M, Adam C, Baulac M, Martinerie J, Varela FJ. 1998. 
% Nonlinear interdependencies of EEG signals in human intracranially 
% recorded temporal lobe seizures. 
% Brain Res. 792:24?40.
% 
% Pereda E, Rial R, Gamundi A, González J. 2001. 
% Assessment of changing interdependencies between human 
% electroencephalograms using nonlinear methods. 
% Phys D Nonlinear Phenom. 148:147?158.


%Inputs:
%x: original time series
%x: predicted time series

Delta = 1 - DPrmse(x,xPred) ./ DPrmse( x , repmat(mean(x,1),1,size(x,1) ) );