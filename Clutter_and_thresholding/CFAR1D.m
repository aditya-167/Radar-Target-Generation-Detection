%T : Number of Training Cells

%G : Number of Guard Cells

%N : Total number of Cells

%Define the number of training cells and guard cells

%Start sliding the window one cell at a time across the complete FFT 1D array. 
%Total window size should be: 2(T+G)+CUT

%For each step, sum the signal (noise) within all the leading or lagging training cells

%Average the sum to determine the noise threshold

%Using an appropriate offset value scale the threshold

%Now, measure the signal in the CUT, which is T+G+1 from the window starting point

%Compare the signal measured in 5 against the threshold measured in 4

%If the level of signal measured in CUT is smaller than the threshold measured, 
%then assign 0 value to the signal within CUT.

clc;
clear all;
close all;

%Data Pts
Ns=1000;

%Generate Random Noise
s=randn(Ns,1);

%Targets location. Assigning bin 100, 200, 300 and 700 as Targets with the amplitudes of 8, 9, 4, 11.
s([100,150,520,600,800,950]) = [ 6, 8, 11, 4, 2, 9];

T_cell = 12;
G_cell = 4;

thresh_cfar = [];
signal_cfar = [];

offset=5;

for i = 1:(Ns-(T_cell+G_cell+1))
    
    noise_level=sum(s(i:i+T_cell-1));
    thresh = (noise_level/T_cell)*offset;
    
    thresh_cfar = [thresh_cfar, {thresh}];
    
    signal=s(i+T_cell+G_cell);
    
    if (signal<thresh)
        signal=0;
    end
    signal_cfar = [signal_cfar, {signal}];
end

plot (cell2mat(signal_cfar),'g--');

% plot original sig, threshold and filtered signal within the same figure.
figure, plot(s);
hold on,plot(cell2mat(circshift(thresh_cfar,G_cell)),'r--','LineWidth',2)
hold on, plot (cell2mat(circshift(signal_cfar,(T_cell+G_cell))),'g--','LineWidth',2);
legend('Signal','CFAR Threshold','detection')
