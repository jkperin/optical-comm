%% Calculate BER for amplified IM-DD link
% bertail = calculate BER using saddlepoint approximation for the tail
% probability. Can't be innacurate because of singularity at origin.
% However, it's much faster
% berpdf (optional) = calculate BER using saddlepoint approximation for the
% pdf, and then calculate tail probability using numerical integration
function [bertail, berpdf] = ber_soa_klse_freq(D, Phi, nu, mpam, tx, soa, rx, sim)

% Maximum frequency (same as used in klse_freq.m)
Fmax = min(1.5*rx.optfilt.fcnorm/2, 0.5);

Mct = round(2*Fmax*sim.Mct); % redefine oversampling ratio to simulate 
Fmax = 1/2*Mct/sim.Mct;
% 'continuous time' so it matches the maximum frequency Fmax
Nsymb = mpam.M^sim.L; % number of data symbols 
Nzero = sim.L; % number of zero symbols pad to begin and end of sequnece.
N = Mct*(Nsymb + 2*Nzero); % total number of points 

% Define time and frequency 
df = 2*Fmax/N;
f = (-Fmax:df:Fmax-df).';

% Generate optical signal
try 
    dataTX = debruijn_generator(mpam.M, sim.L).' - 1;
catch e % can fail for some states of the random number generator
    load rng_state  % state that works
    rng(rng_state);
    dataTX = debruijn_generator(mpam.M, sim.L).' - 1;
    rng('shuffle')
end

% Generate symbols
Pd = mpam.a(gray2bin(dataTX, 'pam', mpam.M) + 1);

% modulate
Pt = reshape(kron(Pd, mpam.pshape(1:Mct)).', [], 1);

% Rescale to desired power
Pmin = 2/(10^(abs(tx.rex)/10)-1); % normalized minimum power when signal has unit average power
% Plevels = (mpam.a/mean(Pt) + Pmin)*soa.Gain*tx.Ptx/(1 + Pmin);
Pthresh = (mpam.b/mean(Pt) + Pmin)*soa.Gain*tx.Ptx/(1 + Pmin);
Pt = (Pt/mean(Pt) + Pmin)*soa.Gain*tx.Ptx/(1 + Pmin); % after rescaling E(Pt) = tx.Ptx

% Calculate electric field (no chirp)
x = sqrt(Pt);
x = [zeros(Nzero*Mct, 1); x; zeros(Nzero*Mct, 1)]; % zero pad
X = fftshift(fft(x));

%%
yy = 0;
xn = zeros(N, sim.Me);
xnd = zeros(Nsymb+2*Nzero, sim.Me);
for k = 1:sim.Me
    phi = Phi(:, k);

    % Interpolate before converting to time domain
    phi = spline(nu, phi, f);
    
    if sim.verbose && k <= 5 
        figure(100), hold on, box on
        plot(f, abs(phi).^2)
        xlabel('f/f_s', 'FontSize', 12)
        ylabel('|\phi_n(f/f_s)|^2', 'FontSize', 12)
    end
    
    % Renormalize
    a = trapz(f, abs(phi).^2);
    phi = 1/sqrt(a)*phi;
    D(k) = sqrt(a)*D(k);
    
    % xn in 'continuous' time
    xn(:, k) = ifft(ifftshift(X.*conj(phi)));
    
    % Xn sampled at decision instants (symbol rate)
    xnd(:, k) = xn(ceil(Mct/2):Mct:end, k);
    
%     yy = yy + D(k).*abs(xn(:, k)).^2;
end

% Remove zeros from begining and end of the sequence
xnd([1:Nzero N/Mct-Nzero+1:N/Mct], :) = [];

%%
varASE = soa.N0*sim.fs/2; % variance of circularly complex Gaussain noise 
varTherm = rx.N0*rx.elefilt.noisebw(sim.fs)/2;

pe = 0;
dat = gray2bin(dataTX, 'pam', mpam.M);
for k = 1:Nsymb
    if dat(k) == mpam.M-1
        pe = pe + tail_saddlepoint_approx(Pthresh(end), D, xnd(k, :).', varASE, varTherm, 'left');
    elseif dat(k) == 0;
        pe = pe + tail_saddlepoint_approx(Pthresh(1), D, xnd(k, :).', varASE, varTherm, 'right');
    else 
        pe = pe + tail_saddlepoint_approx(Pthresh(dat(k) + 1), D, xnd(k, :).', varASE, varTherm, 'right');
        pe = pe + tail_saddlepoint_approx(Pthresh(dat(k)), D, xnd(k, :).', varASE, varTherm, 'left');
    end
end

pe = real(pe)/Nsymb;

bertail = pe/log2(mpam.M);

%% Calculate pdfs using saddlepoint approximation
if nargout == 2
    y = linspace(0, 2*max(abs(x).^2), 1e3);
    px = pdf_saddlepoint_approx(y, D, xnd, varASE, varTherm, sim.verbose);
    if sim.verbose
        for k = 1:length(Pthresh)
            plot(Pthresh(k)*[1 1], [0, 5e4])
        end
    end
    
    pe2 = 0;
    for k = 1:Nsymb
        if dat(k) == mpam.M-1
            pe2 = pe2 + trapz(x(x < Pthresh(end)), px(x < Pthresh(end), k));
        elseif dat(k) == 0
            pe2 = pe2 + trapz(x(x > Pthresh(1)), px(x > Pthresh(1), k));
        else 
            pe2 = pe2 + trapz(x(x > Pthresh(dat(k) + 1)), px(x > Pthresh(dat(k) + 1), k));
            pe2 = pe2 + trapz(x(x < Pthresh(dat(k))), px(x < Pthresh(dat(k)), k));
        end
    end
    pe2 = real(pe2)/Nsymb;
    
    berpdf = pe2/log2(mpam.M);   
end



    