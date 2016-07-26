%% Validate duobinary
clear, clc, close all

addpath ../mpam/
addpath ../f/

Rs = 25e9;
Mct = 10;
Nsymb= 1024;
fs = Mct*Rs;

pulse_shape = select_pulse_shape('rect', 1);

mpam = PAM(4, Rs, 'equally-spaced', pulse_shape);
mpam = mpam.set_levels(0:mpam.M-1, 0.5+(0:mpam.M-2));

dataTX = randi([0 mpam.M-1], [1 Nsymb]);

xd = mpam.signal(dataTX);

xd_enc = duobinary_encoding(xd);

scatterplot(xd_enc-3)

mpam = mpam.norm_levels();

dataRX = mpam.demod(abs(xd_enc-3)/3);

biterr(dataTX, dataRX)

