%% Calculate BER for amplified IM-DD link
function [px, x] = px_soa_klse_freq(Ptx, soa, rx, OptFilt, EleFilt, sim)

sim.Nsymb = 4;
sim.N = sim.Mct*(sim.Nsymb + 2*sim.Nzero);

% Redefine time and frequency 
dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

%
Pt = Ptx*ones(sim.N, 1);

% Calculate electric field (no chirp)
x = sqrt(Pt*soa.Gain);
X = fftshift(fft(x));

% Calculate KL series expansion in the frequency domain
Fmax = 1.5*OptFilt.BW/sim.fs;
[xnd, D] = klse_freq(X, OptFilt.Ho, EleFilt.He, sim, Fmax);

% Calculate pdfs using saddlepoint approximation
x = linspace(0, 5*max(abs(x).^2), 1e3);
px = pdf_saddlepoint_approx(x, D, xnd, soa.N0*sim.fs/2, rx.N0*EleFilt.BW, sim);




    