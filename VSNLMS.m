function [e,w] = VSNLMS(x,d,Lw,mu,psi,alpha,eta,vsFlag)
%VSNLMS Variable Step-size Nomalized LMS algorithm for AEC
%   
w=zeros(Lw,1);
e=d;
Pd=1; Pyhat=1; Pe=1;
for n=length(w):length(x)
    xtdl=x(n:-1:n-Lw+1);
    yhat=w'*xtdl;
    e(n)=d(n)-yhat;
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
    w=w+mu/(xtdl'*xtdl+psi)*xtdl*e(n);
end
end

