function [Pn, CS] = palloc(ofdm, Gch, varNoise, BERtarget)
%% Optimized power allocation and bit loading using Margin-Adaptive Levin-Campello algorithm
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

if strcmpi(class(varNoise), 'function_handle')
    varNoise = varNoise(ofdm.fc);
end

% Gain-to-noise ratio
GNR = abs(Gch).^2./varNoise; 

% Run Levin-Campello algorithm to determine optimal power allocation
[~, CS, Pn] = Levin_Campello_MA(ofdm.B, 1, BERtarget, GNR); 