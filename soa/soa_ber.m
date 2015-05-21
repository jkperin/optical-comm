%% Calulate MGF for the system with SOA
clear, clc, close all

addpath f

% Simulation parameters
sim.Nsymb = 2^14;
sim.Mct = 8;
sim.Nzero = 4;
sim.Me = 32; % Number of eigenvalues
sim.L = 2; % de Bruijin sub-sequence length (ISI symbol length)
sim.M = 10; % Ratio of optical filter BW and electric filter BW (must be integer)
sim.levelSpacing = 'uniform';
sim.verbose = true;
sim.BERtarget = 1e-4;

% M-PAM
mpam.M = 4;
mpam.Rb = 100e9;
mpam.Rs = mpam.Rb/log2(mpam.M);
mpam.pshape = ones(1, sim.Mct);

% Generate unipolar PAM signal
if strcmp(sim.levelSpacing, 'uniform')
    mpam.a = (0:2:2*(mpam.M-1)).';
    mpam.b = (1:2:(2*(mpam.M-1)-1)).';
elseif strcmp(sim.levelSpacing, 'nonuniform')
    %[a, b] = calcOptLevelSpacing(mpam, tx, rx, sim);   
else
    error('Invalide Option!')
end

% 
sim.N = sim.Mct*(sim.Nsymb + 2*sim.Nzero);
sim.fs = mpam.Rs*sim.Mct;

% Time and frequency
dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';
td = t(sim.Mct/2:sim.Mct:end);

sim.t = t;
sim.f = f;

PtxdBm = -20:2:-6;
Ptx = 1e-3*10.^(PtxdBm/10);
tx.rex = 10;  % extinction ratio in dB
tx.Pmin = 1e-6;
rx.N0 = (30e-12).^2;

% Electric Lowpass Filter
EleFilt.type = 'gaussian';
EleFilt.order = 4;
EleFilt.BW = 1.25*mpam.Rs; % single-sided BW
[be, ae] = design_filter(EleFilt.type, EleFilt.order, EleFilt.BW/(sim.fs/2), sim.Mct);
EleFilt.grpdelay = grpdelay(be, ae, 1);
EleFilt.He = @(f) freqz(be, ae, f, sim.fs).*exp(1j*2*pi*f/sim.fs*EleFilt.grpdelay);
% He = @(f) freqz(be, ae, f, sim.fs);

% Optical Bandpass Filter
OptFilt.type = 'gaussian';
OptFilt.order = 4;
OptFilt.BW = sim.M*EleFilt.BW; % single-sided BW
[bo, ao] = design_filter(OptFilt.type, OptFilt.order, OptFilt.BW/(sim.fs/2), sim.Mct);
OptFilt.grpdelay = grpdelay(bo, ao, 1);
OptFilt.Ho = @(f) freqz(bo, ao, f, sim.fs).*exp(1j*2*pi*f/sim.fs*OptFilt.grpdelay);
% Ho = @(f) freqz(bo, ao, f, sim.fs);
 
% SOA
soa = soa(10^(10/10), 20, 1310e-9, 'Chi2', 20);

% [a, b] = soa_level_spacing_optmization(mpam, tx, soa, rx, OptFilt, EleFilt, sim)

for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);
    
    % Mininum power due to finite extinction ratio
    tx.Pmin = 2*tx.Ptx/(10^(tx.rex/10)-1);

    dataTX = randi([0 mpam.M-1], 1, sim.Nsymb);
    
    Pd = mpam.a(gray2bin([zeros(1, sim.Nzero) dataTX zeros(1, sim.Nzero)], 'pam', mpam.M) + 1);
    Pthresh = mpam.b;

    % modulate
    Pt = reshape(kron(Pd, mpam.pshape).', sim.N, 1);
    
    % Scale (!! rescaling should use somethind better than mpam.a)
    Pthresh = (Pthresh/mean(mpam.a) + tx.Pmin/tx.Ptx)*tx.Ptx*soa.Gain;
    Pt = (Pt/mean(mpam.a) + tx.Pmin/tx.Ptx)*tx.Ptx;

    % Calculate electric field (no chirp)
    x = sqrt(Pt);
    
    % Amplifier
    et = soa.amp(x, sim.fs);
    Ef = fftshift(fft(et));

    % Optical bandpass filter
    eo = ifft(fft(et).*ifftshift(OptFilt.Ho(f)));

    % Direct detection and add thermal noise
    yt = abs(eo).^2 + sqrt(rx.N0*sim.fs/2)*randn(size(eo));
    
    % Electric low-pass filter
    yt = real(ifft(fft(yt).*ifftshift(EleFilt.He(f))));

    % Remove zeros from begining and end of the sequence
    nzero = [1:sim.Mct*sim.Nzero sim.N-sim.Mct*sim.Nzero+1:sim.N];
    yt(nzero) = []; 
    Pt(nzero) = []; 

    % Sample
    yd = yt(sim.Mct/2:sim.Mct:end);
    
    % Rescale signal for detection
    yd = (yd - tx.Pmin*soa.Gain)/(tx.Ptx*soa.Gain)*mean(mpam.a);

    % Demodulate
    dataRX = sum(bsxfun(@ge, yd, mpam.b.'), 2);
    dataRX = bin2gray(dataRX, 'pam', mpam.M).';

    % True BER
    [~, ber_mc(k)] = biterr(dataRX, dataTX);
    
    % approx ber
    [ber_approx(k), bertail(k)] = ber_soa_klse_freq(mpam, tx, soa, rx, OptFilt, EleFilt, sim);
    
%     if sim.verbose
%         hold on
%         yd2 = yt(sim.Mct/2:sim.Mct:end);
%         [n, xbin] = hist(yd2(dataTX == 3), 50);
%         n = n/trapz(xbin, n);
%         bar(xbin, n)
%         1;
%     end
%     
end

figure
plot(PtxdBm, log10(ber_mc), PtxdBm, log10(ber_approx))
legend('Monte Carlo', 'Approx')
axis([PtxdBm(1) PtxdBm(end) -10 0])
    