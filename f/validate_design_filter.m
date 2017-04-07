%% Validate design_filter function
clear, clc, close all

fc = 1;
fs = 10;

%% Matlab design
% Butterworth
Bu = design_filter('butter', 5, fc/(fs/2), true);

% Chebychev type 1
C1 = design_filter('cheby1', 5, fc/(fs/2), true);

% Ellyptical
E = design_filter('ellipt', 5, fc/(fs/2), true);

% FIR 1
F1 = design_filter('fir1', 10, fc/(fs/2), true);

%% Impulse invariance transformation of continuous-time transfer function
% Two-pole
P2 = design_filter('two-pole', fc, fs, true);

% Bessel 
B = design_filter('bessel', 5, fc/(fs/2), true);

% Lorentizian
L = design_filter('lorentzian', [], fc/(fs/2), true);

%% FIR approximation
% Gaussian
G1 = design_filter('gaussian', 1, fc/(fs/2), true);
G5 = design_filter('gaussian', 5, fc/(fs/2), true);

% Fiber Bragg gratting
F = design_filter('fbg', [], fc/(fs/2), true);


