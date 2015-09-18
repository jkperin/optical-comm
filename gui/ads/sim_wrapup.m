
dBm2Watt = @(x) 1e-3*10.^(x/10);

load bitfile
save rawsampleData4PAMlong

% 
%% ADS Simulation Parameters
% These values should be provided by ADS simulation
Mct = round(1/(Ts*Rs)); % oversampling ratio
Ndiscard = 100*Mct; % Number of points to discard due to transient

% Interpolate to get odd number of samples per symbol
dt1 = 1/Mct;
dt2 = 1/(Mct + 1);

t1 = 0:dt1:(length(Pout)-1)*dt1;
t2 = 0:dt2:t1(end);

Pout = interp1(t1, Pout, t2);
Poutnf = interp1(t1, Poutnf, t2);
Anode = interp1(t1, Anode, t2);
Mct = Mct + 1;
%
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
ix = Ndiscard+1:length(Pout);
N = length(ix);

xt = xt(ix);
Anode = Anode(ix);
Pout = Pout(ix);
Poutnf = Poutnf(ix);

% Find relative dealy
MaxLag = 20*Mct;
Lags = -MaxLag:MaxLag;
ix = 1:N;
% ix = 1:min(2^12, length(xt)); % enough samples to calculate crosscorrelation
corrAnode = xcorr(xt(ix)-mean(xt(ix)), Anode(ix)-mean(Anode(ix)), MaxLag, 'coeff');
corrPout = xcorr(xt(ix)-mean(xt(ix)), Pout(ix)-mean(Pout(ix)), MaxLag, 'coeff');
corrPoutnf = xcorr(xt(ix)-mean(xt(ix)), Poutnf(ix)-mean(Poutnf(ix)), MaxLag, 'coeff');

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

eyediagram(Pout(1:min(2^14, length(Pout))), 2*Mct)
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

%
dataTXref = mpam.demod(xtref(floor(Mct/2)+1:Mct:end));

save SampleData4PAM

