%% Investigate clipping noise shapping
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
zeroSC = 3; % subcarrier set to zero
% ofdm.CSn(80:end) = 0;

% Generates OFDM symbols
[xncp, symbsTXm] = ofdm.signal(Nsymb);

% Ntaps = 5;
% f = (1:ofdm.Nc/2-1)*ofdm.fc(1)/ofdm.fs;
% Z = zeros(length(f), Ntaps);
% for k = 1:Ntaps
%     Z(:, k) = exp(-1j*2*pi*f.'*k);
% end
% Mag = zeros(size(f));
% Mag(f > ofdm.fc(end)/ofdm.fs) = 1;
% Mag = 1 - Mag.';
% 
% g = pinv(Z)*Mag;

f = 2*(0:ofdm.Nc/2)*ofdm.fc(1)/ofdm.fs;
Mag = zeros(size(f));
Mag(f > 2*ofdm.fc(end)/ofdm.fs) = 1;


% Ntaps = 8;
% bb = fir2(Ntaps, f, Mag);
% % b = bb/bb(1);
% b = bb;
% % b = [1    -1     1    -1     1    -1]; %[1 -1];
% g = -b;
% g(1) = [];
% G = ClassFilter(g, 1, ofdm.fs);
% H = ClassFilter(1, b, ofdm.fs);
% G.removeGroupDelay = false;
% H.removeGroupDelay = false;
% HG = H.cascade(G);
% HG.removeGroupDelay = false;

Ntaps = 5;
[b,a] = yulewalk(Ntaps,f,Mag); % Design IIR filter to meet f, Mag
G = ClassFilter([0 1], 1, ofdm.fs); % G will be delay
H = ClassFilter(1, [1 -1], ofdm.fs); % G will be delay
HG = H; % z^-1 is already implemented within the loop
G.removeGroupDelay = false;
H.removeGroupDelay = false;
HG.removeGroupDelay = false;

sig = sqrt(ofdm.var());

xnf = H.filter(xncp);
H.reset();

yk = zeros(size(xncp));
for t = 1:length(xncp)-1
    
    yk(t+1) = xnf(t) - HG.filter(yk(t));
    
    % Clip
    if abs(yk(t+1)) > rclip*sig
        yk(t+1) = sign(yk(t+1))*rclip*sig;
    end
end

yko = yk;
[yk, lag] = align_signals(yk, xncp);
lag

% Detection
eq.mu = 1e-2;
eq.Ntrain = Inf;
eq.trainSeq = symbsTXm;
[Xn, AGCn, W] = ofdm.detect(yk, eq, true);

% Plots
figure
plot(mean(abs(Xn).^2, 2))

maxLag = 100;
c1 = xcorr(Xn(zeroSC, :), Xn(1, :)-symbsTXm(1, :), maxLag, 'coeff');
c2 = xcorr(Xn(zeroSC, :), Xn(2, :)-symbsTXm(2, :), maxLag, 'coeff');
% c3 = xcorr(Xn(zeroSC, :), Xn(zeroSC, :), maxLag, 'coeff');
c3 = xcorr(Xn(zeroSC, :), Xn(4, :)-symbsTXm(4, :), maxLag, 'coeff');
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