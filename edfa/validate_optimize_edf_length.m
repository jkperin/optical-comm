%% Validate max_channels_on.m
clear, clc, close all

addpath ../f/

E = EDF(10, 'principles_type2');

df = 50e9;
dlamb = df2dlamb(df);
lamb = 1530e-9:dlamb:1565e-9;
SMF = fiber(50e3, @(lamb) 0.18, @(lamb) 0);

Pon = 1e-4;
Signal = Channels(lamb, Pon, 'forward');
Pump = Channels(1480e-9, 50e-3, 'forward');

[~, SpanAttdB] = SMF.link_attenuation(Signal.wavelength);
SpanAttdB = SpanAttdB*ones(size(Signal.wavelength));

[Lopt_interp, Signal_interp] = optimize_edf_length(E, Pump, Signal, Pon, SpanAttdB, 'interp', true);

[Lopt_fminbnd, Signal_fminbnd] = optimize_edf_length(E, Pump, Signal, Pon, SpanAttdB, 'fminbnd', true);
   

