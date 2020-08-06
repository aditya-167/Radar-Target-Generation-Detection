clc;
close all;
clear all;

%R = c*Ts*fb/(2*Bsweep)
%Bsweep = c/2*dr
%Rmax = (Ptransmit*(G^2)*(lambda^2)*(cross-area)/Pmin*4*pi)^1/4
%Tchirp = 5.5*2*Rmax/c

%Note : The sweep time can be computed based on the time needed for the signal to travel 
%the maximum range. In general, for an FMCW radar system, the sweep time should be at least 
%5 to 6 times the round trip time. This example uses a factor of 5.5

Rmax=300;
dr=1;
c=3e8;
fb=[0,1.1e6,13e6,24e6];

Ts=5.5*2*(Rmax/c);

Bsweep=c/(2*dr);
R=c*Ts*fb/(2*Bsweep);

disp(R)