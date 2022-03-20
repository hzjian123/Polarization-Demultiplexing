function [p,t,filtDelay]=raisedCosineFunction(alpha,L,Nsym)%L=8 Naym=8
%Generatea raised cosine (RC) pulse
%alpha - roll-off factor,
%L - oversampling factor
%Nsym - filter span in symbols
%Returns the output pulse p(t) that spans the discrete-time
%base -Nsym:1/L:Nsym. Also returns the filter delay when the
%function is viewed as an FIR filter
Tsym=1; t=-(Nsym/2):1/L:(Nsym/2); % +/- discrete-time base
A = sin(pi*t/Tsym)./(pi*t/Tsym); B=cos(pi*alpha*t/Tsym);
%handle singularities at p(0) and p(t=+/-1/2a)
p = A.*B./(1-(2*alpha*t/Tsym).^2);
p(ceil(length(p)/2))=1; %p(0)=1 and p(0) occurs at center
temp=(alpha/2)*sin(pi/(2*alpha));
p(t==Tsym/(2*alpha))=temp; p(t==-Tsym/(2*alpha))=temp;
filtDelay = (length(p)-1)/2; %FIR filter delay = (N-1)/2
end