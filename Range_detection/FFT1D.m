fs = 1e3;
t=1.5e3
T = 1/fs %sampling period
xtime = (0:t-1)/fs;
F1 = 77;
F2 = 43;
signal=0.7*sin(2*pi*F1*xtime) + 2*sin(2*pi*F2*xtime);
X= signal + 2 *randn(size(xtime));

%plot(1000*xtime(1:80) ,X(1:80))
%title('Signal Corrupted with Zero-Mean Random Noise')
%xlabel('t (milliseconds)')
%ylabel('X(t)')
%%FFT

Y=fft(X);
P2 = abs(Y/t);

%remove mirror image of signal

P1 = P2(1:t/2+1);
f=fs*(0:t/2)/t;

plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
