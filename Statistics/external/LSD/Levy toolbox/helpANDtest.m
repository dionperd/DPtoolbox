% Lévy stable distributions
% 
%    In recent decades, Lévy stable distributions have been playing an increasing role in many diverse scientific and engineering fields. However, it is well known that its analytical expression is largely not available. The following gives a brief introduction to Lévy stable distributions and a new MATLAB toolbox, called LSD, we have developed to evaluate the distributions. You can download this MATLAB toolbox here.
% 
%    The following reference should be cited whenever the ?LSD? is used, ?A survey on computing Lévy stable distributions and a new MATLAB toolbox?, Signal Processing 93 (2013) 242?251. http://dx.doi.org/10.1016/j.sigpro.2012.07.035 (pdf)
% 
% Introduction to Lévy stable distributions
% 
% Lévy stable distributions are generally defined by characteristic functions and a complete specification requires four parameters: stability index, skewness parameter , scale parameter , and location parameter .
% Except of its very few cases such as Gaussian and Cauchy, the closed form expression of density and distribution functions of the Lévy stable distributions is not available.
% The  characteristic function most often employed in numerical calculation is the following:
% 
% Compute pdf and cdf
% 
% There generally are two useful methods, direct integration method and FFT. Direct integration method favors small sample sets since it can be computed at any chosen point, whereas the FFT depends on large data sets to increase accuracy at the expense of higher burdensome computation. Moreover, FFT is only effective for density function with large stability index.
% Generate random numbers
% 
%    Kanter first gave a direct method for Lévy stable random variables, then it was adapted to the general case by Chambers, Mallows and Stuck. Based on the method of Chambers et al., an algorithm of simulating a standard Lévy stable random variable was proposed by Weron , which is regarded as the fastest and the most accurate method.
% 
% Estimate the parameters
% 
% McCulloch obtained four consistent estimators in terms of five sample quantiles and tabulated the values of the four estimators. Empirical characteristic function method is slower, yet much more accurate than quantile estimation.Then, several modifications have been made to this approach, most of all, the implementations of Koutrouvelis , Kogon and Williams  yield good estimators with the restriction. The most accurate but slowest method is ML estimation, which uses either direct integration method or FFT to approximate the density function. However, in many real life problems the higher accuracy does not justify the application of ML estimation, especially when calculations are to be performed on-line.
% MATLAB toolbox
% 
% We choose the relatively reasonable computational algorithms in developing the new MATLAB toolbox to evaluate Lévy stable distributions. When to approximate the density and distribution functions, direct numerical integration introduced by Nolan is utilized efficiently and easy to program. For constructing a standard Lévy stable random variable, the algorithm selected is given by Weron, then can be generalized to all admissible values of the four parameters. Taking the precision and calculation speed into account, we adopt the empirical characteristic function method given by Koutrouvelis to estimate the four parameters of Lévy stable distributions.
% Modeling

%%Example 1: Plot the density function.
x = -10:0.5:10;
semilogy( x , levystblpdf(x, 1.4, 0.8, 1.3, 0.7))
grid on
 

%%Example 2: Plot the distribution function.
x = -10:0.5:10;
semilogy( x , levystblcdf(x, 1.4, 0.8, 1.3, 0.7))
grid on
 
 
%%Example 3: Fit the empirical density of random numbers.
M = levystblrnd(1.4,0.9,1,-0.5,1000);
M=M(:);
x=-20:1:20;
y=hist(M,x);
yy=y/length(M);
bar(x,yy)
hold on
fit = levystblfit(M);
plot(x,levystblpdf(x,fit(1),fit(2),fit(3),fit(4)))
grid on

fit%=
%1.3993 0.8931 0.9999 -0.5058

?
% GUI
% 
% Example:
% Let a random variable , calculate and plot the density and distribution 
% functions in the interval , then generate a random variables matrix with 
% 100 dimensions, and estimate the four parameters.
%                       
% 
% ?
% 
% References:
% 1. Nolan, J. P.: Numerical calculation of stable densities and distribution functions. Commun. Statist. Stochastic Models. 13, 759?774 (1997)
% 2. Nolan, J.P.: An algorithm for evaluating stable densities in Zolotarev?s (M) parameterization. Mathematical and Computer Modeling. 29. 229?233 (1999)
% 3. Weron, R.: On the Chambers-Mallows-Stuck method for simulating skewed stable random variables. Statistics and Probability Letters. 28(2), 165?171 (1996).
% 4. McCulloch, J. H.: Simple consistent estimators of stable distribution parameters. Communications in Statistics-Simulations. 15(4), 1109?1136 (1986)
% 5. Press, S. J.: Estimation in univariate and multivariate stable distribution. Journal of the American Statistical Association. 67(340), 842?846 (1972)
% 6. Koutrouvelis, I.A.: Regression-type estimation of the parameters of stable laws. Journal of the American Statistical Association. 75(372), 918?928 (1980).
% 7. Kogon, S. M., Williams, D. B.: Characteristic function based estimation of stable parameters. in R. Adler, R. Feldman, M. Taqqu (eds.), A Practical Guide to Heavy Tails. Boston: Birkhauser, pp. 311?335 (1998)
% 8. Nolan, J.P.: Stable Distributions: Models for Heavy Tailed Data. Boston: Birkhauser, pp. 3?21 Chapter 1 online at academic2.american.edu/~jpnolan (2010)
% 9. Cízek, P., Weron, R., Härdle, W.: Statistical Tools for Finance and Insurance. Springer?Verlag Berlin Heidelberg. pp. 2?17 (2005)
% 10. Robert, H.R., Nolan, J.P.: Stable Distributions in Mathematica. The Mathematica Journal. 9(4), 776?789 (2005)