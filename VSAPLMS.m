function [e,w] = VSAPLMS(x,d,Lw,mu,alpha,eta,vsFlag)
%VSAPLMS Variable Step-size Affine Projection LMS algorithm for AEC
%   此处显示详细说明
e=d;
M=5;
X=zeros(Lw,M);
w=zeros(Lw,1);
Pd=1; Pyhat=1; Pe=1;
for n=Lw:Lw+M-2
    xtdl=x(n:-1:n-Lw+1);
    X=[xtdl X(:,1:M-1)];
end
for n=Lw+M-1:length(x)-M
    xtdl=x(n:-1:n-Lw+1);
    X=[xtdl X(:,1:M-1)];
    Yhat=X'*w;
    yhat=Yhat(1);
    E=d(n:-1:n-M+1)-Yhat;
    e(n)=E(1);
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
    w=w+mu*X/(X'*X+0.1*eye(M))*E;
end
end

