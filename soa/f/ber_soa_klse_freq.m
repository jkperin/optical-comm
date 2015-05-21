%% Calculate BER for amplified IM-DD link
function [ber, bertail] = ber_soa_klse_freq(mpam, tx, soa,  rx, OptFilt, EleFilt, sim)

sim.L = 2; % de Bruijin sub-sequence length (ISI symbol length)
sim.Nsymb = mpam.M^sim.L;
sim.N = sim.Mct*(sim.Nsymb + 2*sim.Nzero);

% Redefine time and frequency 
dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

% Generate optical signal
rng('shuffle');
try 
    dataTX = debruijn_generator(mpam.M, sim.L).' - 1;
catch e
    load rng_state  % one that works
    rng(rng_state);
    dataTX = debruijn_generator(mpam.M, sim.L).' - 1;
end

Pd = mpam.a(gray2bin([zeros(1, sim.Nzero) dataTX zeros(1, sim.Nzero)], 'pam', mpam.M) + 1);
Pthresh = mpam.b;

% modulate
Pt = reshape(kron(Pd, mpam.pshape).', sim.N, 1);

% Scale (!! rescaling should use somethind better than mpam.a)
Pt = (Pt/mean(mpam.a) + tx.Pmin/tx.Ptx)*tx.Ptx;

% Calculate electric field (no chirp)
x = sqrt(Pt*soa.Gain);
X = fftshift(fft(x));

Pthresh = (Pthresh/mean(mpam.a) + tx.Pmin/tx.Ptx)*tx.Ptx*soa.Gain;

% Calculate KL series expansion in the frequency domain
Fmax = 1.5*OptFilt.BW/sim.fs;
[xnd, D] = klse_freq(X, OptFilt.Ho, EleFilt.He, sim, Fmax);

% Calculate pdfs using saddlepoint approximation
x = linspace(0, 2*max(abs(x).^2), 1e3);
px = pdf_saddlepoint_approx(x, D, xnd, soa.N0*sim.fs/2, rx.N0*EleFilt.BW, sim);
if sim.verbose
    for k = 1:length(Pthresh)
        plot(Pthresh(k)*[1 1], [0, 5e4])
    end
end

pe = 0;
pe2 = 0;
dat = gray2bin(dataTX, 'pam', mpam.M);
for k = 1:sim.Nsymb
    if dat(k) == mpam.M-1
%         pe = pe + tail_saddlepoint_approx(Pthresh(end), D, xnd(k, :).', soa.N0*sim.fs/2, rx.N0*EleFilt.BW, 'left');
        pe2 = pe2 + trapz(x(x < Pthresh(end)), px(x < Pthresh(end), k));
        
    elseif gray2bin(dataTX(k), 'pam', mpam.M) == 0;
%         pe = pe + tail_saddlepoint_approx(Pthresh(1), D, xnd(k, :).', soa.N0*sim.fs/2, rx.N0*EleFilt.BW, 'right');
        pe2 = pe2 + trapz(x(x > Pthresh(1)), px(x > Pthresh(1), k));
        
    else 
%         pe = pe + tail_saddlepoint_approx(Pthresh(dat(k) + 1), D, xnd(k, :).', soa.N0*sim.fs/2, rx.N0*EleFilt.BW, 'right');
%         pe = pe + tail_saddlepoint_approx(Pthresh(dat(k)), D, xnd(k, :).', soa.N0*sim.fs/2, rx.N0*EleFilt.BW, 'left');
        pe2 = pe2 + trapz(x(x > Pthresh(dat(k) + 1)), px(x > Pthresh(dat(k) + 1), k));
        pe2 = pe2 + trapz(x(x < Pthresh(dat(k))), px(x < Pthresh(dat(k)), k));
    end
end


pe = real(pe)/sim.Nsymb;
pe2 = real(pe2)/sim.Nsymb;

% Use tail probability based on pdf and integration because it is more
% accurate than tail probabibility based on saddle-point approximation
ber = pe2/log2(mpam.M);
bertail = pe/log2(mpam.M);



    