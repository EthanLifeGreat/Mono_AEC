function [e,w] = VSNLMSNt(x,d,Lw,mu,mu_po,epsilon,psi,gamma,alpha,beta,eta,vsFlag)
%VSNLMSNT Variable Step-size Nomalized LMS-Newton algorithm for AEC
%   此处显示详细说明
N=Lw; M=8;
P=ones(M+1,1);
w=zeros(N,1);
kappa=zeros(M,1);
b=x(N+M:-1:N);
bp=zeros(M,1);
u_a=zeros(N,1);
f=zeros(M,1);
fp=zeros(M,1);
e=d;
Pd=1; Pyhat=1; Pe=1;
for n=N+M+1:length(x)
    %
    %	Lattice Predictor
    %
    b_old=b;
    f(1)=x(n);
    b(1)=f(1);
    P(1)=beta*P(1)+0.5*(1-beta)*(f(1)^2+b_old(1)^2);
    for m=1:M
        f(m+1)=f(m)-kappa(m)*b_old(m);
        b(m+1)=b_old(m)-kappa(m)*f(m);
        kappa(m)=kappa(m)+2*mu_po*(P(m)+epsilon)^(-1)*(f(m)*b(m+1)+b_old(m)*f(m+1));
        P(m+1)=beta*P(m+1)+0.5*(1-beta)*(f(m+1)^2+b_old(m+1)^2);
        if abs(kappa(m))>gamma
            kappa(m)=sign(kappa(m))*gamma;
        end
    end
    %
    %  u_a(n) update
    %
    u_a(2:N)=u_a(1:N-1);
    b_old=bp;
    fp(1)=b(M+1);
    bp(1)=fp(1);
    for m=1:M
        fp(m+1)=fp(m)-kappa(m)*b_old(m);
        bp(m+1)=b_old(m)-kappa(m)*fp(m);
    end
    u_a(1)=fp(M+1)/(P(M+1)+epsilon);
    %
    %   Filtering
    %
    xtdl=x(n-M:-1:n-N-M+1);
    yhat=w'*xtdl;
    e(n)=d(n-M)-yhat;
    % Step-size Variation
    if vsFlag==1
        Pd=alpha*Pd+(1-alpha)*d(n)*d(n);
        Pyhat=alpha*Pyhat+(1-alpha)*yhat*yhat;
        Pe=alpha*Pe+(1-alpha)*e(n)*e(n);
        mu=1-eta*Pyhat/Pd;
        if mu>1
            mu=1;
        elseif mu<0
            mu=0;
        end
    end
    if (u_a'*xtdl)>=0
        w=w+(mu/((u_a'*xtdl)+psi))*e(n)*u_a;
    end
end
end