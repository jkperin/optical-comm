%% Discrete-time simulation of Stokes vector receiver
clear, clc, close all

addpath f/
addpath ../f/
addpath ../mpam/
addpath ../apd/

sim.Nsymb = 2^15;
sim.Rb = 111e9;
rexdB = -15;
Prx = 1e-3;
BERtarget = 1.8e-4;
Nrealizations = 1e4;
amplified = false;

%% Modulation format specified for each polarization
% PAM(M, bit rate, leve spacing : {'equally-spaced', 'optimized'}, pulse
% shape: struct containing properties of pulse shape 
pulse_shape = select_pulse_shape('rect', 1);
IM = PAM(4, sim.Rb/3, 'equally-spaced', pulse_shape);
PM = PAM(4, sim.Rb/3, 'equally-spaced', pulse_shape);
PM.a = [-1/2; 0; 1/2; 1];
PM.b = [-1/4; 1/4; 3/4];
sim.fs = IM.Rs;

%% ========================== Amplifier ===================================
% Constructor: OpticalAmplifier(Operation, param, Fn, Wavelength)
% - Opertation: either 'ConstantOutputPower' or 'ConstantGain'
% - param: GaindB if Operation = 'ConstantGain', or outputPower
% if Operation = 'ConstantOutputPower'
% - Fn:  noise figure in dB
% - Wavelength: operationl wavelength in m
OptAmp = OpticalAmplifier('ConstantOutputPower', 0, 5, 1550e-9); 
% OptAmp = OpticalAmplifier('ConstantGain', 20, 5, Tx.Laser.wavelength);
OptAmp.maxGaindB = 50;

%% ============================ Receiver ==================================
%% Photodiodes
% pin(R: Responsivity (A/W), Id: Dark current (A), BW: bandwidth (Hz))
% PIN frequency response is modeled as a first-order system
Rx.PD = pin(1, 10e-9, Inf);

%% TIA-AGC
% One-sided thermal noise PSD
Rx.N0 = (30e-12).^2; 

%% ============================ Hybrid ====================================
% polarization splitting --------------------------------------------------
Rx.PolSplit.sig  = 'PBS';                                                   % pbs: polarization beamsplitter
Rx.PolSplit.Rext = 30;                                                      % PBS extinction ratio (dB), default = 30

%% Transmission
PrxdBmreq = zeros(1, Nrealizations);
OSNRdBreq = zeros(1, Nrealizations);
dataTX = randi([0 IM.M-1], [3, sim.Nsymb]); % symbol stream for each polarization
AB = zeros(2, Nrealizations);
for n = 1:Nrealizations
    phi = rand(1, 3)*2*pi;
    U1 = [exp(-1j*phi(1)/2), 0; 0 exp(1j*phi(1)/2)];
    U2 = [cos(phi(2)/2) -1j*sin(phi(2)/2); -1j*sin(phi(2)/2) cos(phi(2)/2)];
    U3 = [cos(phi(3)/2) -sin(phi(3)/2); sin(phi(3)/2) cos(phi(3)/2)];
    U = U1*U2*U3; % = [a -b; b^* a^*], for a and b complex

    a = U(1, 1);
    b = -U(1, 2);
    
    AB(1, n) = a;
    AB(2, n) = b;

    M = [abs(a)^2, abs(b)^2, -2*real(a*b'), 2*imag(a*b');...
        abs(b)^2, abs(a)^2,  2*real(a*b'), -2*imag(a*b');...
        real(a*b), -real(a*b),  real(a^2) - real(b^2), -imag(a^2) - imag(b^2);...
        imag(a*b), -imag(a*b), imag(a^2) - imag(b^2), real(a^2) + real(b^2)];

    Minv = inv(M);

%     U = eye(2);
    if amplified
        PtxdBm = -35:-20;
    else
        PtxdBm = -15:0;
    end
    PrxdBm = zeros(size(PtxdBm));
    BER = zeros(3, length(PtxdBm));
    for k = 1:length(PtxdBm)
        IM = IM.adjust_levels(dBm2Watt(PtxdBm(k)-3), rexdB);

        Px = IM.mod(dataTX(1, :));
        Py = IM.mod(dataTX(2, :));
        phixy = exp(1j*pi*PM.mod(dataTX(3, :)));

        Etx = [sqrt(Px).*phixy; sqrt(Py)];
        
        Erx = U*Etx;
        if amplified
            [Erec, OSNRdB(k)] = OptAmp.amp(Erx, sim.fs);
        else
            Erec = Erx;
        end
        
        PrxdBm(k) = power_meter(Erec);
        Y = dual_pol_stokes_receiver(Erec, Rx, sim);
        Y = 2*Y;
        Y = Minv*Y;

        phiXY = angle(Y(3, :) + 1j*Y(4, :))/pi; % {-pi, -pi/2, 0, pi/2, pi} -> {-1, -1/2, 0, 1/2, 1}
        phiXY(phiXY < -3/4) = 1; % resolve -pi, pi ambiguity
        
        IM = IM.adjust_levels(dBm2Watt(PrxdBm(k)-3), rexdB);
        dataRX(1, :) = IM.demod(Y(1, :));
        dataRX(2, :) = IM.demod(Y(2, :));
        dataRX(3, :) = PM.demod(phiXY);
        
        % Measure noise
%         [~, Ps] = power_meter(Erec);
%         W1 = Y(1, :) - IM.mod(dataTX(1, :));
%         histfit(W1)
%         var(W1)
%         2*Ps*OptAmp.Ssp*IM.Rs/2
%         
%         W2 = Y(2, :) - IM.mod(dataTX(2, :));
%         histfit(W2)
%         var(W2)
%         2*Ps*OptAmp.Ssp*IM.Rs/2
%         
%         W3 = Y(3, :) - real(sqrt(IM.mod(dataTX(1, :))).*exp(1j*pi*PM.mod(dataTX(3, :))).*sqrt(IM.mod(dataTX(2, :))));
%         histfit(W3)
%         var(W3)
%         Ps*OptAmp.Ssp*IM.Rs/2
%         
%         W4 = Y(4, :) - imag(sqrt(IM.mod(dataTX(1, :))).*exp(1j*pi*PM.mod(dataTX(3, :))).*sqrt(IM.mod(dataTX(2, :))));
%         histfit(W4)
%         var(W4)
%         Ps*OptAmp.Ssp*IM.Rs/2
        
        [~, BER(1, k)] = biterr(dataTX(1, :), dataRX(1, :));
        [~, BER(2, k)] = biterr(dataTX(2, :), dataRX(2, :));
        [~, BER(3, k)] = biterr(dataTX(3, :), dataRX(3, :));
    end
    
%     figure(100), hold on, box on
%     plot(PtxdBm, log10(BER))
%     plot(PtxdBm, log10(mean(BER, 1)))
%     legend('BER X', 'BER Y', 'BER ang(XY)', 'BER')
    
    BERtotal = mean(BER, 1);
    if amplified
        OSNRdBreq(n) = fit_ber(OSNRdB, BERtotal, BERtarget);
    else
        PrxdBmreq(n) = fit_ber(PrxdBm, BERtotal, BERtarget);
    end
        
end

figure
if amplified
    boxplot(OSNRdBreq)
else
    boxplot(PrxdBmreq)
end
    