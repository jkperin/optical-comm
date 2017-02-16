%% Investigate correlation of clipping noise across subcarriers
% Conclusions:
% Clipping noise is approximately Gaussian distributed in each subcarrier
% The clipping noise in a given subcarrier is uncorrelated to the clipping
% noise in neighboring subcarriers
% Clipping noise in a given subcarrier is apparently white 
clear, clc, close all

addpath ../
addpath ../f/
addpath ../../f/

Nsymb = 2^10;
rclip = 1;

% Create OFDM class
ofdm = ofdm(256, 208, 16, 112e9);
ofdm.set_cyclic_prefix(0, 0);
ofdm.CSn = ofdm.CS*ones(1, ofdm.Nu/2);
zeroSC = 21; % subcarrier set to zero
ofdm.CSn(zeroSC:2:30) = 16;

% Generates OFDM symbols
[xncp, symbsTXm] = ofdm.signal(Nsymb);

sig = sqrt(ofdm.var());

% Clip
figure, hold on, box on
plot(xncp)
plot([1 length(xncp)], rclip*sig*[1 1], ':k')
plot([1 length(xncp)], -rclip*sig*[1 1], ':k')

yk = xncp;
yk(xncp > rclip*sig) = rclip*sig;
yk(xncp < -rclip*sig) = -rclip*sig;

% Detection
[Xn, AGCn, W] = ofdm.detect(yk, [], true);

% Plots
figure
plot(mean(abs(Xn).^2, 2))

scatterplot(Xn(zeroSC, :))

maxLag = 100;
c1 = xcorr(Xn(zeroSC, :), Xn(1, :)-symbsTXm(1, :), maxLag, 'coeff');
c2 = xcorr(Xn(zeroSC, :), Xn(2, :)-symbsTXm(2, :), maxLag, 'coeff');
% c3 = xcorr(Xn(zeroSC, :), Xn(zeroSC, :), maxLag, 'coeff');
c3 = xcorr(Xn(zeroSC, :), Xn(5, :)-symbsTXm(5, :), maxLag, 'coeff');
figure, hold on
plot(-maxLag:maxLag, abs(c1).^2)
plot(-maxLag:maxLag, abs(c2).^2)
plot(-maxLag:maxLag, abs(c3).^2)
title('Correlation coefficient of clipping noise in neighboring subcarriers')

%% Clipping noise power
clipNoiseE2 = mean(abs(Xn-symbsTXm).^2, 2);
clipNoiseVar = var(Xn-symbsTXm, 0, 2);
figure, box on, hold on
stem(1:ofdm.Nu/2, clipNoiseE2)
stem(1:ofdm.Nu/2, clipNoiseVar)
xlabel('Subcarrier index')
ylabel('Clipping noise variance')
legend('Second moment', 'variance')

    
