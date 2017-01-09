function xP = DPpredAvTransl(u,Np,nnInd,x,X,k,D,ind_i,ind_i_plus_u)


%This function calculates the average predicted translation

% Schiff SJ, So P, Chang T, Burke RE, Sauer T. 1996. 
% Detecting dynamical interdependence and generalized synchrony through 
% mutual prediction in a neural ensemble. 
% Phys Rev E Stat Phys Plasmas Fluids Relat Interdiscip Topics.
% 54:6708?6724.

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
%u: prediction horizon in sample points
%Np: number of predicted points
%nnInd: indexes of k nearest neighbors of each point in the X pointset in a
%   matrix of (points,neighbors)
%x: univariate signal
%X: embedded pointset of x
%k: number of nearest neighbors
%D: dimensionality of the embedding
%ind_i: indices of points that predict
%ind_i_plus_u: indices of points to be predicted

%Outputs:
%xP: the predicted points


clear X k D;

%Preallocate memory
xP = zeros(Np,1);

%Predict each point...
for iP=1:Np;
    %...as the average translation of each nearest neighbors
    xP(iP) = mean( x( nnInd(iP,:) + u ) );
end