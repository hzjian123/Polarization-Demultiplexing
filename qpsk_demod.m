function [a_cap,x,y] = qpsk_demod(r,fc,OF)
%Function to demodulate a conventional QPSK signal
%r - received signal at the receiver front end
%fc - carrier frequency in Hertz
%OF - oversampling factor (multiples of fc) - at least 4 is better
%L - upsampling factor on the inphase and quadrature arms
%a_cap - detected binary stream
fs = OF*fc; %sampling frequency
L = 2*OF; %samples in 2Tb duration
t=0:1/fs:(length(r)-1)/fs; %time base
x=r.*cos(2*pi*fc*t); %I arm
y=-r.*sin(2*pi*fc*t); %Q arm
x = conv(x,ones(1,L));%integrate for L (Tsym=2*Tb) duration
y = conv(y,ones(1,L));%integrate for L (Tsym=2*Tb) duration
x = x(L:L:end);%I arm - sample at every symbol instant Tsym
y = y(L:L:end);%Q arm - sample at every symbol instant Tsym
x=x/OF;y=y/OF;
a_cap = zeros(1,2*length(x));
a_cap(1:2:end) = x.' ;%> 0; %even bits
a_cap(2:2:end) = y.' ;%> 0; %odd bits
 %To plot constellation at the receiver
