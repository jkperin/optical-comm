%% Delay vs bandwidth tradeoff (theory)
clear, clc, close all

addpath ../f/
addpath ../../f/

Rb = 2*112e9;
sim.ModFormat = QAM(4, Rb/2);
M = 4;
Rs = Rb/(2*log2(M));
BERtarget = 1.8e-4;
sim.BERtarget = BERtarget;
SNRdBRef = SNRreq(BERtarget, M, 'QAM');
SNRRef = 10^(SNRdBRef/10);
linewidth = 2*200e3;
csi = 1/sqrt(2);
Ts = 1/sim.ModFormat.Rs;
varPN = 2*pi*linewidth*Ts;
etac = mean(abs(sim.ModFormat.mod(0:M-1)).^2)/mean(1./abs(sim.ModFormat.mod(0:M-1)).^2);

wn = 2*pi*linspace(0, 300, 100)*1e6;
w = 2*pi*linspace(-0.5*sim.ModFormat.Rs, 0.5*sim.ModFormat.Rs, 1e3);
Delayps = (0:50:800)*1e-12;

figure(1), hold on, box on
for k = 1:length(Delayps)
    GammaPN = zeros(size(wn));
    GammaAWGN = zeros(size(wn));
    phiError = zeros(size(wn));
       
    Delay = Delayps(k);
    for kk = 1:length(wn)
        % Calculate phase error
        numFs = [2*csi*wn(kk) wn(kk)^2];
        denFs = [1 0];

        Fw = polyval(numFs, 1j*w)./polyval(denFs, 1j*w);

        % 1/(s + Fs) = denFs/(denFs*s + numFs)
        %     [num1, den1] = tfdata(1/(s + Fs));
        H1 = 1j*w + exp(-1j*w*Delay).*Fw;
        GammaPN(kk) = 2*csi*wn(kk)/pi*trapz(w, 1./abs(H1).^2);

        % Fs/(s + Fs) = numFs*denFs/(denFs*s + numFs)
        %     [num2, den2] = tfdata(Fs/(s + Fs));
        H2 = Fw./(1j*w + exp(-1j*w*Delay).*Fw);
        GammaAWGN(kk) = 2*csi/(pi*(1 + 4*csi^2)*wn(kk))*trapz(w, abs(H2).^2);

        % Phase error
        phiError(kk) = varPN/(2*csi*wn(kk)*Ts)*GammaPN(kk) + (1 + 4*csi^2)*wn(kk)*Ts/(4*csi)*(etac/(2*SNRRef))*GammaAWGN(kk);

        Perror = ber_qpsk_imperfect_cpr(SNRdBRef, phiError(kk)^2);
        
        if isinf(Perror) || isnan(Perror)
            SNRdB(k, kk) = NaN;
        else
            SNRdB(k, kk) = SNRreq(Perror, M, 'QAM');
        end
    end
    
    figure(1)
    plot(wn/(2*pi*1e6), SNRdB(k, :), 'DisplayName', sprintf('Delay = %.1f ps', Delayps(k)*1e12))
    drawnow
    
    figure(2), hold on
    plot(wn/(2*pi*1e6), phiError, 'DisplayName', sprintf('Delay = %.1f ps', Delayps(k)*1e12))
    
end

figure(1), legend('-DynamicLegend')
axis([0 wn(end)/(2*pi*1e6) 0 15])
figure(2), legend('-DynamicLegend')
axis([0 wn(end)/(2*pi*1e6) 0 0.03])


