%% Raised cosine pulse shape
clear, clc, close all

addpath ../f/

Nsymb = 1024;
Mct = 11;
M = 4;


%% Pulse shape
pulse_shape.type = 'rc'; % either 'rect', 'rrc', or 'rc'
pulse_shape.sps = Mct; % number of samples per symbol
switch lower(pulse_shape.type) 
    case 'rect' % Rectangular pulse shape
        pulse_shape.h = ones(1, pulse_shape.sps); % pulse shapping filter coefficients
    case 'rrc' % Root raised cosine pulse shape
        pulse_shape.rolloff = 0.25; % 
        pulse_shape.span = 6; % number of symbols over which pulse shape spans
        pulse_shape.h = rcosdesign(pulse_shape.rolloff, pulse_shape.span, pulse_shape.sps, 'sqrt'); % pulse shapping filter coefficients
    case 'rc' % Raised cosine pulse shape
        pulse_shape.rolloff = 0.25; % 
        pulse_shape.span = 6; % number of symbols over which pulse shape spans
        pulse_shape.h = rcosdesign(pulse_shape.rolloff, pulse_shape.span, pulse_shape.sps, 'normal'); % pulse shapping filter coefficients
    otherwise
        error('Invalid pulse shape type')
end


dataTX = randi([0 M-1], [1 Nsymb]);
x = upsample(dataTX, Mct);
xt = filter(pulse_shape.h, 1, x);

figure
n = length(pulse_shape.h);
stem(-(n-1)/2:(n-1)/2, pulse_shape.h)

figure, clf
eyediagram(xt, 2*Mct)