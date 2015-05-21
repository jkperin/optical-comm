%% Calulate MGF for the system with SOA
clear, clc, close all

addpath f

% M-PAM
mpam.M = 4;
mpam.Rb = 100e9;
mpam.Rs = mpam.Rb/log2(mpam.M);

% Simulation parameters
sim.Me = 32; 
sim.Nzero = 4;
sim.L = 2; % de Bruijin sub-sequence length (ISI symbol length)
sim.M = 2; % Ratio of optical filter BW and electric filter BW (must be integer)
sim.Mct = 8;
sim.verbose = true;

sim.Nsymb = mpam.M^sim.L;
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

tx.rex = 10;  % extinction ratio in dB
rx.N0 = (30e-12).^2;

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
 
% SOA
soa = soa(5, 20, 1310e-9, 'Chi2', 20);

% Generate optical signal
Ptx = 1e-5;
Pmin = 2*Ptx/(10^(tx.rex/10)-1);
rng('shuffle');
try 
    dataTX = debruijn_generator(mpam.M, sim.L).' - 1;
catch e
    rng('shuffle');
    dataTX = debruijn_generator(mpam.M, sim.L).' - 1;
end

Pd = pammod([zeros(1, sim.Nzero) dataTX zeros(1, sim.Nzero)], mpam.M, 0, 'gray');
Pthresh = 1:2:mpam.M+1;
pt = ones(1, sim.Mct); % pulse shape
Pt = reshape(kron(Pd, pt).', sim.N, 1);
Pt = Pt + (mpam.M - 1);
Pthresh = (Pthresh/mean(Pt) + Pmin/Ptx)*Ptx;
Pt = (Pt/mean(Pt) + Pmin/Ptx)*Ptx*soa.Gain;
x = sqrt(Pt);
X = fftshift(fft(x));

% Noise
w = sqrt(soa.N0*sim.fs/2)*(randn(sim.N, 1) + 1j*randn(sim.N, 1));

et = x + w;
Ef = fftshift(fft(et));

eo = ifft(fft(et).*ifftshift(Ho(f)));

yt = real(ifft(fft(abs(eo).^2).*ifftshift(He(f))));

% Remove zeros from begining and end of the sequence
nzero = [1:sim.Mct*sim.Nzero sim.N-sim.Mct*sim.Nzero+1:sim.N];
yt(nzero) = []; 
Pt(nzero) = []; 

Fmax = 1.5*OptFilt.BW/sim.fs;

[xnd, D, yyd] = klse_freq(X, Ho, He, sim, Fmax);

x = linspace(0, 2*max(Pt), 100);
px = pdf_saddlepoint_approx(x, D, xnd, soa.N0*sim.fs/2, rx.N0*EleFilt.BW, sim);
for k = 1:length(Pthresh)
    plot(Pthresh(k)*[1 1], [0, 5e4])
end

% for k = 1:length(x);
%     px2(k) = tail_saddlepoint_approx(x(k), D, xnd(1, :).', soa.N0*sim.fs/2, rx.N0*EleFilt.BW, 'right');
% end
% 
% figure
% plot(x, px2, x, cumtrapz(x, px))

pe = 0;
pe2 = 0;
for k = 1:sim.Nsymb
    if dataTX(k) == mpam.M-1
        pe = pe + tail_saddlepoint_approx(Pthresh(end), D, xnd(k, :).', soa.N0*sim.fs/2, rx.N0*EleFilt.BW, 'left');
        pe2 = pe2 + trapz(x(x < Pthresh(end)), px(x < Pthresh(end), k));
        
    elseif dataTX(k) == 0;
        pe = pe + tail_saddlepoint_approx(Pthresh(1), D, xnd(k, :).', soa.N0*sim.fs/2, rx.N0*EleFilt.BW, 'right');
        pe2 = pe2 + trapz(x(x > Pthresh(1)), px(x > Pthresh(1), k));
        
    else 
        pe = pe + tail_saddlepoint_approx(Pthresh(dataTX(k) + 1), D, xnd(k, :).', soa.N0*sim.fs/2, rx.N0*EleFilt.BW, 'right');
        pe = pe + tail_saddlepoint_approx(Pthresh(dataTX(k)), D, xnd(k, :).', soa.N0*sim.fs/2, rx.N0*EleFilt.BW, 'left');
        pe2 = pe2 + trapz(x(x > Pthresh(dataTX(k) + 1)), px(x > Pthresh(dataTX(k) + 1), k));
        pe2 = pe2 + trapz(x(x < Pthresh(dataTX(k))), px(x < Pthresh(dataTX(k)), k));
    end
end

pe = real(pe)/sim.Nsymb;
pe2 = real(pe2)/sim.Nsymb;

%% Figures
figure, hold on
tnz = t(sim.Mct*sim.Nzero+1:end-sim.Mct*sim.Nzero);
plot(tnz, yt)
plot(tnz, Pt)
stem(td(sim.Nzero+1:end-sim.Nzero), yyd, '*')



    