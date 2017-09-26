%% Compute the capacity of Stokes vector receiver when the input is uniformly distributed
clear, clc, close all

addpath ../f/

N = 2^26;
SNRdB = 0:20;
SNR = 10.^(SNRdB/10);

Ndim = 4;
MQAM = QAM(16, 1);

Ex = MQAM.mod(randi([0 MQAM.M-1], [1 N])).';
Ey = MQAM.mod(randi([0 MQAM.M-1], [1 N])).';
Etx = [Ex Ey];
X = [real(Ex) imag(Ex) real(Ey) imag(Ey)];
Ps = MQAM.Pqam;

Nrealizations = 1;
MI = zeros(Nrealizations, length(SNR));
for n = 1:Nrealizations
    % Generate random unitary matrix
    phi = rand(1, 3)*2*pi;
    U1 = [exp(-1j*phi(1)/2), 0; 0 exp(1j*phi(1)/2)];
    U2 = [cos(phi(2)/2) -1j*sin(phi(2)/2); -1j*sin(phi(2)/2) cos(phi(2)/2)];
    U3 = [cos(phi(3)/2) -sin(phi(3)/2); sin(phi(3)/2) cos(phi(3)/2)];

    U = U1*U2*U3; % = [a -b; b^* a^*], for a and b complex
    
    Erx = U*Etx.';
%     Xrx = [real(Erx(1, :)); imag(Erx(1, :)); real(Erx(2, :)); imag(Erx(2, :))]; 
    Xrx = [abs(Erx(1, :)).^2;
        abs(Erx(2, :)).^2;
        2*real(Erx(1, :).*conj(Erx(2, :)));...
        2*imag(Erx(1, :).*conj(Erx(2, :)))];
    for k = 1:length(SNR) 
        Ps = mean(abs(Xrx(:)).^2);
        Yrx = Xrx + sqrt(Ps/(SNR(k)))*randn(4, N);
              
%         Y = U\Yrx;
        Y = Yrx;
        MI(n, k) = mutual_info(X(:, 1:Ndim), Y(1:Ndim, :).', 30);
    end
end

figure(1), box on, hold on
plot(SNRdB, MI, 'LineWidth', 2, 'DisplayName', sprintf('%d-QAM coherent', MQAM.M))
plot(SNRdB, log2(1 + SNR), 'LineWidth', 2, 'DisplayName', 'Shannon 2 DOF')
plot(SNRdB, 2*log2(1 + SNR), 'LineWidth', 2, 'DisplayName', 'Shannon 4 DOF')
xlabel('SNR (dB)')
ylabel('Spectral efficiency (bits/s/Hz)')
legend('-dynamiclegend')
grid on


SNRe = linspace(SNR(1), SNR(end));
SNRdBe = 10*log10(SNRe(1:end-1));
dSNR = abs(SNRe(2)-SNRe(1));

EDOF = @(SNR, MI) SNR(1:end-1).*log(2).*diff(MI)*abs(SNR(2)-SNR(1));

MIfit = spline(SNR, MI, SNRe);

figure(2), box on, hold on
plot(SNRdBe, EDOF(SNRe, MIfit), 'LineWidth', 2, 'DisplayName', sprintf('%d-QAM coherent', MQAM.M))
plot(SNRdBe, EDOF(SNRe, log2(1 + SNRe)), 'LineWidth', 2, 'DisplayName', 'Shannon 2 DOF')
plot(SNRdBe, EDOF(SNRe, 3/2*log2(1 + SNRe)), 'LineWidth', 2, 'DisplayName', 'Shannon 3 DOF')
xlabel('SNR (dB)')
ylabel('EDOF')
legend('-dynamiclegend')
grid on