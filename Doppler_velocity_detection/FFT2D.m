clear all;
close all;
clc;

P=peaks(20)
M=1024
N=128
signal = repmat(P,[M N]);
signal_fft = fft2(signal,1024,128);
signal_fft = fftshift(signal_fft);
signal_fft = abs(signal_fft);
imagesc(signal_fft);
