function [e,wF] = VSNPFBLMS(x,d,Lw,M,mu,psi,alpha,eta,vsFlag)
%VSNPFBLMS Variable Step-size Nomalized Partitioned Frequency-domain/Fast...
%   Block LMS algorithm for AEC
%   
N=Lw; P=N/M;
wF=zeros(2*M,P);
xF=zeros(2*M,P);
e=d;
Pd=1; Pyhat=1; Pe=1;
for n=M+1:M:length(x)-M
    xF=[fft(x(n-M:n+M-1)) xF(:,(1:end-1))];
    yhat=ifft(sum((wF.*xF).').'); 
    yhat=real(yhat(M+1:end));
    E=d(n:n+M-1)-yhat;
    % Step-size Variation
    if vsFlag==1
        Pd=alpha*Pd+(1-alpha)*sum(abs(d(n:n+M-1)));
        Pyhat=alpha*Pyhat+(1-alpha)*sum(abs(yhat));
        Pe=alpha*Pe+(1-alpha)*sum(abs(e(n:n+M-1)));
        mu=1-eta*Pyhat/Pd;
        if mu>1
            mu=1;
        elseif mu<0
            mu=0;
        end
    end
    MU=mu.*(sum((abs(xF).^2),2)+psi).^(-1);
    EF=fft([zeros(M,1); E]);
    wF=wF+diag(MU.*EF)*conj(xF);
    waux=real(ifft(wF)); wF=fft([waux(1:M,:); zeros(M,P)]);
    e(n:n+M-1)=d(n:n+M-1)-yhat;
end

end

