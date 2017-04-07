%% Investigate correlation of clipping noise across subcarriers
% Conclusions:
% Clipping noise is approximately Gaussian distributed in each subcarrier
% The clipping noise in a given subcarrier is uncorrelated to the clipping
% noise in neighboring subcarriers
% Clipping noise in a given subcarrier is apparently white 
clear, close all

addpath ../
addpath ../f/
addpath ../../f/

Nsymb = 2^10;
rclip = 2;
K = 1 - 2*qfunc(rclip);

% Create OFDM class
ofdm = ofdm(256, 208, 16, 112e9);
ofdm.set_cyclic_prefix(0, 0);
ofdm.CSn = ofdm.CS*ones(1, ofdm.Nu/2);
zeroSC = 21; % subcarrier set to zero
% ofdm.CSn(zeroSC:31) = 0;
% ofdm.Pn(zeroSC:31) = 2;

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

d1 = yk - K*xncp;
xk2 = xncp - d1/K;

yk2 = xk2;
yk2(xk2 > rclip*sig) = rclip*sig;
yk2(xk2 < -rclip*sig) = -rclip*sig;

sigx2 = sqrt(sig^2 + var(d1)/K^2);
r2 = rclip*sig/sigx2; 
K2 = 1 - 2*qfunc(r2);

d2 = yk2 - K2*xk2;
xk3 = xk2 - d2/K2;

yk3 = xk3;
yk3(xk3 > rclip*sig) = rclip*sig;
yk3(xk3 < -rclip*sig) = -rclip*sig;

% Detection
[Xn, AGCn, W] = ofdm.detect(yk2, [], true);
ber = ofdm.count_ber([1:128 Nsymb-128:Nsymb])

% maxLag = 100;
% c1 = xcorr(Xn(zeroSC, :), Xn(1, :)-symbsTXm(1, :), maxLag, 'coeff');
% c2 = xcorr(Xn(zeroSC, :), Xn(2, :)-symbsTXm(2, :), maxLag, 'coeff');
% % c3 = xcorr(Xn(zeroSC, :), Xn(zeroSC, :), maxLag, 'coeff');
% c3 = xcorr(Xn(zeroSC, :), Xn(5, :)-symbsTXm(5, :), maxLag, 'coeff');
% figure, hold on
% plot(-maxLag:maxLag, abs(c1).^2)
% plot(-maxLag:maxLag, abs(c2).^2)
% plot(-maxLag:maxLag, abs(c3).^2)
% title('Correlation coefficient of clipping noise in neighboring subcarriers')

%% Clipping noise power
clipNoiseE2 = mean(abs(Xn-K*symbsTXm).^2, 2);
clipNoiseVar = var(Xn-K*symbsTXm, 0, 2);
% figure, box on, hold on
% stem(1:ofdm.Nu/2, clipNoiseE2)
% stem(1:ofdm.Nu/2, clipNoiseVar)
% xlabel('Subcarrier index')
% ylabel('Clipping noise variance')
% legend('Second moment', 'variance')


figure, box on, hold on
stem(1:ofdm.Nu/2, 10*log10(var(Xn, 0, 2)./clipNoiseVar-1))
xlabel('Subcarrier index')
ylabel('Signal-to-clipping noise ratio (dB)')
    
