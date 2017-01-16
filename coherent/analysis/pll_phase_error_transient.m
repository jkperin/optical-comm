%% Transient phase error response in 2nd-order PLL
% This analysis doesn't include the possibility of cycle slips. Hence, the
% actual transient time is much longer than the one shown here
% This analysis also doesn't include loop delay

clear, clc, close all

addpath ../../f/

Rb = 2*112e9;
Npol = 2;
M = 4;
Rs = Rb/(Npol*log2(M));
wn = 2*pi*50e6;
xi = 1/sqrt(2);
N = 2^12;
fs = 50e9;
[f, t] = freq_time(N, fs);

foff = [0.1e9 1e9 2e9 3e9 4e9 5e9];

figure, hold on, box on
leg = {};
for k = 1:length(foff)
    plot(t*Rs, 2*pi*foff(k)/wn*(1/sqrt(1 - xi^2))*sin(sqrt(1-xi^2))*wn*t.*exp(-xi*wn*t))
    leg = [leg ['f_{off} = ' sprintf('%.2f GHz', foff(k)/1e9)]];
end

xlabel('Symbols')
ylabel('Phase error (rad)')
legend(leg)
title('Phase error transient for frequency step in 2nd-order PLL')

