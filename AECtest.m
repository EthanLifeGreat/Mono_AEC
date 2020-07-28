%% AECtest
%   This is a demo for calling various AEC
%	(Acoustic Echo Cancellation) functions, namely
%     - VSNLMS:
%         Variable Step-size Nomalized LMS algorithm
%     - VSNLMSNt: 
%         Variable Step-size Nomalized LMS-Newton algorithm
%     - VSAPLMS:
%         Variable Step-size Affine Projection LMS algorithm
%     - VSNPFBLMS: 
%         Variable Step-size Nomalized Partitioned Frequency-domain/Fast
%         Block LMS algorithm
%     - SbLMS:
%         Subband LMS algorithm
%
%	These AEC functions are copied and modified from 
%	Behrouz Farhang-Boroujeny (2013) Adaptive Filters:Theory and Applications 
%
%	Author: Ethan Chen
%	Last updated: 7/2020

clear; clc; close all;clear sound;

%% Options for Adaptive Filter
% Algorithm Option
%   1.VSNLMS    2.VSNLMSNt  3.VSAPLMS 
%   4.VSNPFBLMS             5.SbLMS
opt.AEC = 5;

% Length of adaptive filter
opt.Lw=1024;

%% Create Input Signals x & d
% Load far-end signal x (reference signal) and it's Smapling rate
load speech3
% Load near-end simulated echo filter w0
load w0
% Length of w0
Lw0=length(w0);
% Time index
t=(0:length(x)-1)/Fs;
% Convolution of x and w to create echo signal y
y=filter(w0,1,x);
% Gaussian random noise e0
e0=.001*randn(size(y));
% Desired signal d
d=y+e0;

%% Call AEC Function
tic
if opt.AEC == 1             % VSNLMS
    mu=0.5; psi=0.1; alpha=0.995; vsFlag=1; eta=0.5;
    [e,w] = VSNLMS(x,d,opt.Lw,mu,psi,alpha,eta,vsFlag);
    
elseif opt.AEC == 2         % VSNLMSNt
    mu=0.5; psi=0.1; alpha=0.995; vsFlag=1; eta=0.5;
    mu_po=0.001; epsilon=0.05; beta=0.98; gamma=0.9;
    [e,w] = VSNLMSNt(x,d,opt.Lw,mu,mu_po,epsilon,psi,gamma,alpha,beta,eta,vsFlag);
    
elseif opt.AEC == 3         % VSAPLMS
    % NOTE: This is EXTREMELY time-consuming.
    mu=0.5; psi=0.1; alpha=0.995; vsFlag=1; eta=0.5;
    [e,w] = VSAPLMS(x,d,opt.Lw,mu,alpha,eta,vsFlag);
    
elseif opt.AEC == 4         % VSNPFBLMS
    mu=0.5; psi=0.1; alpha=0.995; vsFlag=1; eta=0.5;
    M=64;
    [e,wF] = VSNPFBLMS(x,d,opt.Lw,M,mu,psi,alpha,eta,vsFlag);
    
elseif opt.AEC == 5         % SbLMS
    L=16;
    [e,w,h_a,h_s] = SbLMS(x,d,opt.Lw,L);
    
else 
    error(['Invalid opt.AEC value: ', int2str(opt.AEC)])   
end
toc

%% Check Output
%sound(e,Fs);
if opt.AEC == 1 || opt.AEC == 2 || opt.AEC == 3 % VSNLMS&VSNLMSNt&VSAPLMS
    xhat = filter(1,w,y);
    yhat = filter(w,1,x);
    figure(2)
    subplot(421),plot(w0),title("echo filter");
    subplot(422),plot(w),title("adaptive filter");
    subplot(423),plot(t,x),title("reference signal");
    subplot(424),plot(t,d),title("desired signal");
    subplot(425),plot(t,xhat),title("recovered signal");
    subplot(426),plot(t,yhat),title("simulated desired signal");
    subplot(4,2,[7,8]),plot(t,d,'c'),hold on,plot(t,e,'b'),hold off,title("residual during learning");
    
elseif opt.AEC == 4  % VSNPFBLMS
    yhat = PFBfilter(wF,x);
    figure(2)
    subplot(421),plot(w0),title("echo filter");
    %
    subplot(423),plot(t,x),title("reference signal");
    subplot(424),plot(t,d),title("desired signal");
    % 
    subplot(426),plot(t,yhat),title("simulated desired signal");
    subplot(4,2,[7,8]),plot(t,d,'c'),hold on,plot(t,e,'b'),hold off,title("residual during learning");
    
else  % SbLMS
    yhat = SbFilter(w,x,L,h_a,h_s);
    figure(2)
    subplot(421),plot(w0),title("echo filter");
    %
    subplot(423),plot(t,x),title("reference signal");
    subplot(424),plot(t,d),title("desired signal");
    %
    subplot(426),plot(t,yhat),title("simulated desired signal");
    subplot(4,2,[7,8]),plot(t,d,'c'),hold on,plot(t,e,'b'),hold off,title("residual(IFFTed from frequency domain)");
    
end