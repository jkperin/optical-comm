%% Optical power penalty due to extinction ratio as a function of average optical power

clear, clc, close all

% DC-OFDM
pp_dc = @(rex) 10*log10(10^(rex/10)+1) - 10*log10(10^(rex/10)-1);

% ACO-OFDM
pp_aco = @(rex, r) 10*log10(10^(rex/10)+(r*sqrt(2*pi) - 1)) - 10*log10(10^(rex/10)-1);

rex = 10; % dB
fprintf('rex = %d\nDC-OFDM: %.5f\nACO-OFDM: %.5f\n', rex, pp_dc(rex), pp_aco(rex+3,4))

rex = 20; % dB
fprintf('rex = %d\nDC-OFDM: %.5f\nACO-OFDM: %.5f\n', rex, pp_dc(rex), pp_aco(rex+3,4))
