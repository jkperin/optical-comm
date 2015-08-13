clear, clc, close all

CodesPath = 'C:\Users\jose.krauseperin\Documents\codes';

addpath([CodesPath '\f\'])
addpath([CodesPath '\mpam\'])

dBm2Watt = @(x) 1e-3*10.^(x/10);

load bitfile

% 
%% ADS Simulation Parameters
% These values should be provided by ADS simulation
% M = 4;
% Ts = 2.5e-12;   % Sampling period in ADS simulation
% Rs = 25e9;      % Symbol rate

Mct = round(1/(Ts*Rs)); % oversampling ratio
Ndiscard = 20*Mct; % Number of points to discard due to transient

if M == 2
    symbolStream = bitStream;
elseif M == 4
    % ADS simulation must be set to little Endian with Gray encoding
    tmp = reshape(bitStream, log2(M), []);
    symbolStream = tmp(2, :)*2 + tmp(1, :); % converts to decimal
%     symbolStream = bin2gray(symbolStream, 'pam', M); % converts to gray
else
    error('PAM order not supported.')
end

% Generate ideal PAM signal
pshape = @(n) double(n >= 0 & n < Mct); % pulse shape
mpam = PAM(M, Rs*log2(M), 'equally-spaced', pshape);
xt = mpam.mod(symbolStream, Mct);

% Discard samples compromised in the transient
N = length(Pout)-Ndiscard;
ix = Ndiscard+1:N;

xt = xt(ix);
Anode = Anode(ix);
Pout = Pout(ix);
Poutnf = Poutnf(ix);

% Find relative dealy
MaxLag = 20*Mct;
Lags = -MaxLag:MaxLag;
corrAnode = xcorr(xt-mean(xt), Anode-mean(Anode), MaxLag, 'coeff');
corrPout = xcorr(xt-mean(xt), Pout-mean(Pout), MaxLag, 'coeff');
corrPoutnf = xcorr(xt-mean(xt), Poutnf-mean(Poutnf), MaxLag, 'coeff');

figure, hold on, grid on, box on
plot(Lags, corrAnode)
plot(Lags, corrPout)
plot(Lags, corrPoutnf)
xlabel('Lag')
ylabel('Correlation Coefficient')
title('Correlation between data sequence and ADS signals')
legend('Anode', 'Pout', 'Pout w/o filter')

% Align symbol sequence and ADS signals
[~, l1] = max(abs(corrAnode));
l1 = Lags(l1);
Anode = circshift(Anode, [0 l1]);

[~, l2] = max(abs(corrPout));
l2 = Lags(l2);
Pout = circshift(Pout, [0 l2]);

[~, l3] = max(abs(corrPoutnf));
l3 = Lags(l3);
Poutnf = circshift(Poutnf, [0 l3]);

% Discard samples from begining and end
l = max(abs(l1), abs(l2));
Anode([1:l*Mct end-l*Mct+1:end]) = [];
Pout([1:l*Mct end-l*Mct+1:end]) = [];
Poutnf([1:l*Mct end-l*Mct+1:end]) = [];
xt([1:l*Mct end-l*Mct+1:end]) = [];

eyediagram(Pout, 2*Mct)
Pthresh = [4.9000    7.7000   10.5000];
hold on
for k = 1:length(Pthresh)
    plot([-0.5 0.5], Pthresh(k)*[1 1], 'k')
end

% Get integer number of symbols
Nsymb = floor(length(Anode)/Mct);
N = Nsymb*Mct;
Anode = Anode(1:N);
Pout = Pout(1:N);
Poutnf = Poutnf(1:N);
xtref = xt(1:N);

save SampleData4PAMlong

