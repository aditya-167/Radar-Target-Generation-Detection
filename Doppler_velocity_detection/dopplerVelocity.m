clc;
close all;
clear all;

%Fd: shift in transmitted frequency due to doppler
%vr = velcity of target
%lambda = wavelength of signal
%Fd = 2*vr/lambda

%The beat frequency comprises of both frequency components: f_rf 
%(frequency delta due to range) and f_df 
%(frequency shift due to velocity). Although, in the case of automotive radar the f_df 
%is very small in comparison to the f_rf 

%Hence, the doppler velocity is calculated by measuring the rate of change of phase across multiple chirps.
%Δφ= delx/lambda
%Δφ= f*delx/c
%del(f)=del(phi)/del(t)

frequency=77e9;
c=3e8;
lambda=c/frequency;
doppler_f = [3e3,4.5e3,11e3,-3e3];

velocity = doppler_f*lambda/2;
disp(velocity)
