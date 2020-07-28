function [e,w,h_a,h_s,yhat] = SbLMS(x,d,Lw,L)
%SBLMS Subband LMS algorithm for AEC
%   
% parameters
K_a=5;
K_s=3;
M_a=19;
M_s=32;M=M_s;
Delta_a=K_a*M_a;
Delta_s=K_s*M_s;
N_a=2*K_a*M_a+1;
N_s=2*K_s*M_s+1;
%Delta=ceil((Delta_a+Delta_s)/M)*M;
alpha_a=3/16;
alpha_s=7/19;

%%
% filter design
h_a=Eigenfir(M_a,N_a,alpha_a,Delta_a)';
h_s=Eigenfir(M_s,N_s,alpha_s,Delta_s)';

Lw_sub=ceil(Lw/L);
w=zeros(Lw_sub,M/2);
e=0;
yhat=0;
mu=0.5;
%%
for i=1:M/2
    h_a_i=h_a.*exp(1j*(2*pi*(i-1)/M)*(0:length(h_a)-1));
    h_s_i=h_s.*exp(1j*(2*pi*(i-1)/M)*(0:length(h_s)-1));
    x_sub=decimate(filter(h_a_i,1,x),L);
    d_sub=decimate(filter(h_a_i,1,d),L);
    w_sub=zeros(Lw_sub,1);
    e_sub=zeros(size(d_sub));
    y_sub_hat=zeros(size(d_sub));
    for n=Lw_sub:length(x_sub)
        xtdl=x_sub(n:-1:n-Lw_sub+1);
        y_sub_hat(n)=w_sub'*xtdl;
        e_sub(n)=d_sub(n)-y_sub_hat(n);
        w_sub=w_sub+mu/(xtdl'*xtdl+0.00001)*xtdl*e_sub(n)';
    end
    e=e+real(filter(h_s_i,1,expander(e_sub,L)));
    yhat=yhat+real(filter(h_s_i,1,expander(y_sub_hat,L)));
    w(:,i)=w_sub;
end
yhat=2*L*real(yhat);
yhat=yhat(1:length(x));
e=e(1:length(x));
end

%
% EXPANDER: y=expander(x,L) 
% When x is a vector, this function adds L-1 after each element of x. 
% When x is a matrix, each column of it treated as vector and expanded L
% fold.
function y=expander(x,L)
[M,N]=size(x);
if (N==1)||(M==1)
    if M<N
        y=zeros(1,N*L);
        y(1:L:end)=x;
    else
        y=zeros(M*L,1);
        y(1:L:end)=x;
    end
else
    y=zeros(M*L,N);
    y(1:L:end,:)=x;
end
end

%*******************************************************
%** This function can be added to MATLAB's vocabulary **
%** for designing the propotype filter of an M-band   **
%** complementary filter bank with the controlable    **
%** delay. The input parameters are:                  **
%**    M     ---  No. of bands of the filter bank;    **
%**    N     ---  Length of the designed filter;      **
%**    alpha ---  Roll-off factor;                    **
%**    del   ---  Group delay, i.e., Delta.           **
%** The output vector h stores the coefficients of    **
%** the designed FIR filter.                          **
%*******************************************************         
%
%   Format:       h = eigenfir(M,N,alpha,del)
% 
%
% Last updated on April 28, 1998
%

function h = Eigenfir(M,N,alpha,del)
  h=zeros(N,1);
  
  ws = (1+alpha)/M;  
  corr_row = (-1)*ws*sinc(ws*(0:N-1));
  corr_row(1) = 1-ws;
  corr = toeplitz(corr_row);

  left0 = (del-M:-M:1); left0 = left0(size(left0,2):-1:1); 
  right0 = (del+M:M:N); 
  zer0 = [left0 right0];
  corr(zer0,:) = []; corr(:,zer0) = [];

  [vect,eigen] = eig(corr); 
  [~,mcnu] = min(diag(real(eigen)));
  vh = vect(:,mcnu); 
  j = 1; k = 1; 
  
  for i = 1:N
     if (i == zer0(j)) 
       h(i) = 0; 
       if j<size(zer0,2)
	 j = j+1;
       end 
     else
       h(i) = vh(k); 
       k = k+1; 
     end
  end 
  h = h/h(del)/M;
end