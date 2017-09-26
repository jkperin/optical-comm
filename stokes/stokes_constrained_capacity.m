%% Compute the capacity of Stokes vector receiver when the input is uniformly distributed
clear, clc

addpath ../mpam/
addpath ../f/

N = 2^26;
SNRdB = -5:30;
SNR = 10.^(SNRdB/10);

Ndim = 4;
IM = PAM(8, 1);
IM.a = IM.a + 0.1;
PM = PAM(8, 1);
% PM.a = [-1/2; 0; 1/2; 1];
d = 2/PM.M;
PM.a = (-1+d:d:1).';

Px = IM.mod(randi([0 IM.M-1], [1 N]));
Py = IM.mod(randi([0 IM.M-1], [1 N]));
phixy = exp(1j*pi*PM.mod(randi([0 PM.M-1], [1 N])));
Etx = [sqrt(Px).*phixy; sqrt(Py)];

clear Px Py phixy % clear variable to free space

Pspam = IM.Ppam;

Nrealizations = 1;
MI = zeros(Nrealizations, length(SNR));
tic
for n = 1:Nrealizations
    % Generate random unitary matrix
    phi = rand(1, 3)*2*pi;
    U1 = [exp(-1j*phi(1)/2), 0; 0 exp(1j*phi(1)/2)];
    U2 = [cos(phi(2)/2) -1j*sin(phi(2)/2); -1j*sin(phi(2)/2) cos(phi(2)/2)];
    U3 = [cos(phi(3)/2) -sin(phi(3)/2); sin(phi(3)/2) cos(phi(3)/2)];

    U = U1*U2*U3; % = [a -b; b^* a^*], for a and b complex

    a = U(1, 1);
    b = -U(1, 2);

    % Define fiber transfer matrix in Muller space
    % As defined in [1] M. Chagnon et al,"Digital signal processing for
    % dual-polarization intensity and interpolarization phase modulation 
    % formats using stokes detection", J. Lightw. Technol. 34, 188–195. 2016.
    M = [abs(a)^2, abs(b)^2, -2*real(a*b'), 2*imag(a*b');...
        abs(b)^2, abs(a)^2,  2*real(a*b'), -2*imag(a*b');...
        real(a*b), -real(a*b),  real(a^2) - real(b^2), -imag(a^2) - imag(b^2);...
        imag(a*b), -imag(a*b), imag(a^2) - imag(b^2), real(a^2) + real(b^2)];

    X = [abs(Etx(1, :)).^2;
        abs(Etx(2, :)).^2;
        2*real(Etx(1, :).*conj(Etx(2, :)));...
        2*imag(Etx(1, :).*conj(Etx(2, :)))].';
    
    Xrx = M*X.';
       
    for k = 1:length(SNR) 
        Ps = mean(abs(Xrx(:)).^2);
        Y = M\(Xrx + sqrt(Ps/(SNR(k)))*randn(4, N));
%         Y = Xrx + sqrt(Ps/(SNR(k)))*randn(4, N);
        MI(n, k) = mutual_info(X(:, 1:Ndim), Y(1:Ndim, :).', 30);
    end
end
toc

figure(1), box on, hold on
plot(SNRdB, MI, 'LineWidth', 2, 'DisplayName', sprintf('2 x %d-PAM', IM.M))
plot(SNRdB, 0.5*log2(1 + SNR), 'LineWidth', 2, 'DisplayName', 'Shannon 1 DOF')
plot(SNRdB, log2(1 + SNR), 'LineWidth', 2, 'DisplayName', 'Shannon 2 DOF')
plot(SNRdB, 3/2*log2(1 + SNR), 'LineWidth', 2, 'DisplayName', 'Shannon 3 DOF')
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
plot(SNRdBe, EDOF(SNRe, MIfit), 'LineWidth', 2, 'DisplayName', sprintf('2 x %d-PAM', IM.M))
plot(SNRdBe, EDOF(SNRe, log2(1 + SNRe)), 'LineWidth', 2, 'DisplayName', 'Shannon 2 DOF')
plot(SNRdBe, EDOF(SNRe, 3/2*log2(1 + SNRe)), 'LineWidth', 2, 'DisplayName', 'Shannon 3 DOF')
xlabel('SNR (dB)')
ylabel('EDOF')
legend('-dynamiclegend')
grid on


clear Etx X Xrx Y 
save stokes_3x8.mat