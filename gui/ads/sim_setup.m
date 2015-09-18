%% Matlab Setup
function [out1, out2] = sim_setup(in1, in2, in3)

persistent n

if isempty(n)
    n = 1;
end

global Anode Pout Poutnf

% CodesPath = 'C:\Users\jose.krauseperin\Documents\codes';
% 
% addpath([CodesPath '\f\'])
% addpath([CodesPath '\mpam\'])
% addpath([CodesPath '\ofdm\'])
% addpath([CodesPath '\apd\'])
% addpath([CodesPath '\soa\'])

% Mct = 16;
% 
% mpam = PAM(2, 50e9, 'equally-spaced', @(n) double(n >= 0 & n < Mct)); 
% 
% mpammod = @(x) 3*mpam.mod(x, Mct) - 1.5;
% 
% out = [];

out1 = in1;
out2 = in2;

Pout(n) = in2;
Poutnf(n) = in1;
Anode(n) = in3;

n = n + 1;
