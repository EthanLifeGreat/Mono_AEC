function yhat = PFBfilter(wF,x)
%PFBFILTER Uses original signal and output wF of VSNPFBLMS()...
%   to estimate the desired signal
%   
    M=size(wF,1)/2;
    P=size(wF,2);
    xF=zeros(2*M,P);
    yhat=zeros(size(x));
    for n=M+1:M:length(x)-M
        xF=[fft(x(n-M:n+M-1)) xF(:,(1:end-1))];
        yHat=ifft(sum((wF.*xF).').'); 
        yHat=real(yHat(M+1:end));
        yhat(n:n+M-1)=yhat(n:n+M-1)+yHat;
    end
end

