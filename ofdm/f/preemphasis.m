function [Pn, CS] = preemphasis(ofdm, Gch, varNoise, BERtarget)
%% Power allocation using preemphasis method: inverse of channel response and uniform bit loading
% Quantization noise is assumed negligible compared to other noises
% Inputs:
% - ofdm: OFDM class
% - Gch: frequency response of the channel at the subcarriers frequency
% - varNoise: function handle or vector containing noise variance at each
% subacarrier
% - BERtarget: target BER
% Outputs:
% - Pn: power allocation at each subcarrier
% - CS: constellation size at each subcarrier

% Definitions
Nu = ofdm.Nu;
CS = ofdm.CS*ones(1,Nu/2);

if strcmpi(class(varNoise), 'function_handle')
    varNoise = varNoise(ofdm.fc);
end

% SNR to achieve target BER sim.BERtarget
snrdB = fzero(@(x) berqam(ofdm.CS, x) - BERtarget, 20);
snrl = 10^(snrdB/10);

Pnrx = snrl*(varNoise);                                               

Pn = Pnrx./(abs(ofdm.K*Gch).^2);       

% Discard carriers which requrie infinite energy
CS(isinf(Pn)) = 0;
Pn(isinf(Pn)) = 0;
