%% Calulate MGF for the system with SOA
clear, clc, close all

addpath f

% M-PAM
mpam.M = 4;
mpam.Rb = 100e9;
mpam.Rs = mpam.Rb/log2(mpam.M);
% sim.Me = 64;

% Simulation parameters
sim.M = 1; % Ratio of optical filter BW and electric filter BW (must be integer)
sim.Mct = 4;
sim.N = 2^16;
sim.fs = mpam.Rs*sim.Mct;

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

% Electric Lowpass Filter
EleFilt.type = 'gaussian';
EleFilt.order = 4;
EleFilt.BW = 1.25*mpam.Rs; % single-sided BW
[be, ae] = design_filter(EleFilt.type, EleFilt.order, EleFilt.BW/(sim.fs/2), sim.Mct);
EleFilt.grpdelay = grpdelay(be, ae, 1);
He = @(f) freqz(be, ae, f, sim.fs).*exp(1j*2*pi*f/sim.fs*EleFilt.grpdelay);
% He = @(f) freqz(be, ae, f, sim.fs);

% Optical Bandpass Filter
OptFilt.type = 'gaussian';
OptFilt.order = 4;
OptFilt.BW = sim.M*EleFilt.BW; % single-sided BW
[bo, ao] = design_filter(OptFilt.type, OptFilt.order, OptFilt.BW/(sim.fs/2), sim.Mct);
OptFilt.grpdelay = grpdelay(bo, ao, 1);
Ho = @(f) freqz(bo, ao, f, sim.fs).*exp(1j*2*pi*f/sim.fs*OptFilt.grpdelay);
% Ho = @(f) freqz(bo, ao, f, sim.fs);
% 
tx.PtxdBm = -30:-5;
tx.Ptx = 1e-3*10.^(tx.PtxdBm/10);
tx.lamb = 1310e-9;
soa.Fn = 9; % dB
soa.Seq = @(G) (G - 1)*10^(soa.Fn/10)/2*(1.98644568e-25/tx.lamb); % one-sided PSD
soa.G = 10^(5/10);

% Generate optical signal
Nzero = 32;
Ndis = Nzero/sim.Mct;
Nsymb = sim.N/sim.Mct;
dataTX = randi([0 mpam.M-1], [Nsymb 1]);
Pd = pammod(dataTX, mpam.M, 0, 'gray');
pt = ones(1, sim.Mct);
P = reshape(kron(Pd, pt).', sim.N, 1);
P = P + (mpam.M - 1);
Pbar = mpam.M - 1;

for k = 1:length(tx.Ptx)
    
    Pt = P/Pbar*(tx.Ptx(k)*soa.G);
    Pt(1:Nzero) = 0;
    Pt(end-Nzero+1:end) = 0;
    
    x = sqrt(Pt);

    % Noise
    w = sqrt(soa.Seq(soa.G)*sim.fs/2)*(randn(sim.N, 1) + 1j*randn(sim.N, 1));

    et = x + w;

    eo = ifft(fft(et).*ifftshift(Ho(f)));

    yt = real(ifft(fft(abs(eo).^2).*ifftshift(He(f))));
    
    % Sampling
    yd = yt(sim.Mct/2:sim.Mct:end);
    
    yd = yd/(tx.Ptx(k)*soa.G)*Pbar - (mpam.M - 1);
    
    dataRX = pamdemod(yd, mpam.M, 0, 'gray');
    
    [~, ber(k)] = biterr(dataRX(Ndis+1:end-Ndis-1), dataTX(Ndis+1:end-Ndis-1));
end
    
figure
plot(tx.PtxdBm, log10(ber))


% Filters 
% figure
% plot(f/sim.fs, abs(He(f)).^2, f/sim.fs, abs(Ho(f)).^2)
% legend('electric', 'optical')