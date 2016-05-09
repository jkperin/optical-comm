function [Eout1, Eout2] = powerSplitter(Ein1, Ein2, splittingRatio, phaseShiftDeg, timeDelay, f)
%% Optical power splitter
% Example of transfer matrix if 3dB power splitter: [Eout1; Eout2] = 1/sqrt(2)*[1 j; j 1][Ein1; Ein2]
% Inputs:
% - Ein1, Ein2: electric field in input 1 and 2, respectively
% - splittingRatio: power splitting ratio: sqrt(splittingRatio) from port 1
% goes to output 1, and 1-sqrt(splittingRatio) goes to output 2. For port 2
% the power is reversed
% - phaseShiftDeg: additional phase shift between input 1 and input 2. Ein2
% will be delay by phaseShiftDeg in the output
% - timeDelay (optional): time delay in seconds between input 1 and input
% 2. Ein2 will be delayed by timeDelay
% - f (optinal): frequency vector. Only necessary if timeDelay ~= 0

if nargin == 4 || timeDelay == 0
    phaseShift = deg2rad(phaseShiftDeg);
    Eout1 = sqrt(splittingRatio)*Ein1 + sqrt(1-splittingRatio)*1j*Ein2*exp(1j*phaseShift);
    Eout2 = sqrt(1-splittingRatio)*1j*Ein1 + sqrt(splittingRatio)*Ein2*exp(1j*phaseShift);
else
    f = ifftshift(f);
    if size(Ein1, 1) == 2 % two pols
        f = [f; f];
    end
    phaseShift = deg2rad(phaseShiftDeg) - (2*pi*f)*timeDelay;
    Eout1 = sqrt(splittingRatio)*Ein1 + sqrt(1-splittingRatio)*1j*ifft(fft(Ein2, length(Ein2), 2).*exp(1j*phaseShift), length(Ein2), 2);
    Eout2 = sqrt(1-splittingRatio)*1j*Ein1 + sqrt(splittingRatio)*ifft(fft(Ein2, length(Ein2), 2).*exp(1j*phaseShift), length(Ein2), 2);
end