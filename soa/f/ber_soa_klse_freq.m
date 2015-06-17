%% Calculate BER for amplified IM-DD link
% bertail = calculate BER using saddlepoint approximation for the tail
% probability. Can't be innacurate because of singularity at origin.
% However, it's much faster
% berpdf (optional) = calculate BER using saddlepoint approximation for the
% pdf, and then calculate tail probability using numerical integration
function [bertail, berpdf] = ber_soa_klse_freq(D, Phi, Fmax, mpam, tx, soa, rx, sim)

% Maximum frequency (same as used in klse_freq.m)
[nu, ~] = lgwt(sim.Me, -Fmax, Fmax);
nu = nu(end:-1:1); % put in right order

% 'continuous time' so it matches the maximum frequency Fmax
Nsymb = mpam.M^sim.L; % number of data symbols 
Nzero = sim.L; % number of zero symbols pad to begin and end of sequnece.
N = sim.Mct*(Nsymb + 2*Nzero); % total number of points 

% Define time and frequency 
df = 1/N;
f = (-0.5:df:0.5-df).';
fm = f(abs(f) <= Fmax);

% Generate optical signal
dataTX = debruijn_sequence(mpam.M, sim.L).';

% Generate symbols
Pd = mpam.a(gray2bin(dataTX, 'pam', mpam.M) + 1);

% modulate
Pt = reshape(kron(Pd, mpam.pshape(1:sim.Mct)).', [], 1);

% Rescale to desired power
if strcmp(mpam.level_spacing, 'uniform') % adjust levels to include extinction ratio penalty
    Pmin = mean(Pt)/(10^(abs(tx.rex)/10)-1); % minimum power 
    Plevels = (mpam.a + Pmin)*soa.Gain*tx.Ptx/(mean(Pt) + Pmin); % levels at the receiver
    Pthresh = (mpam.b + Pmin)*soa.Gain*tx.Ptx/(mean(Pt) + Pmin); % decision thresholds at the receiver
    Pt = (Pt + Pmin)*soa.Gain*tx.Ptx/(mean(Pt) + Pmin); % after rescaling E(Pt) = tx.Ptx
elseif strcmp(mpam.level_spacing, 'nonuniform') % already includes extinction ratio penalty, so just scale
    Plevels = mpam.a*soa.Gain*tx.Ptx/mean(Pt);
    Pthresh = mpam.b*soa.Gain*tx.Ptx/mean(Pt);
    Pt = Pt*soa.Gain*tx.Ptx/mean(Pt);
end 

% Calculate electric field (no chirp)
x = sqrt(Pt);
x = [zeros(Nzero*sim.Mct, 1); x; zeros(Nzero*sim.Mct, 1)]; % zero pad
X = fftshift(fft(x));

%%
yk = 0;
xn = zeros(N, sim.Me);
for k = 1:sim.Me
    phi = Phi(:, k);

    % Interpolate before converting to time domain
    phi = interp1(nu, phi, fm, 'pchip');
    
    % Renormalize
    a = trapz(fm, abs(phi).^2);
    phi = 1/sqrt(a)*phi;
    
    if sim.verbose && k <= 5
        figure(101)
        subplot(211), hold on, box on
        plot(fm, abs(phi).^2)
        xlabel('f/f_s', 'FontSize', 12)
        ylabel('|\phi_n(f/f_s)|^2', 'FontSize', 12)
        
        subplot(212), hold on, box on
        plot(fm, unwrap(angle(phi)))
        xlabel('f/f_s', 'FontSize', 12)
        ylabel('phase', 'FontSize', 12)
    end
    
    phi = [zeros(sum(f < -Fmax), 1); phi; zeros(sum(f > Fmax), 1)];
    
    % xn in 'continuous' time
    xn(:, k) = ifft(ifftshift(X.*conj(phi)));
       
    yk = yk + D(k).*abs(xn(:, k)).^2;
end
% Discard zeros and downsample
ix = sim.Mct*Nzero+sim.Mct/2:sim.Mct:N-sim.Mct*Nzero;
yd = yk(ix);
xnd = xn(ix, :);

if sim.verbose
    figure(102), hold on
    plot([zeros(Nzero*sim.Mct, 1); Pt; zeros(Nzero*sim.Mct, 1)])
    plot(yk, '-k')
    plot(ix, yd, 'o')
    legend('Transmitted power', 'KL-SE Fourier')
end

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
    y = linspace(0, 2*max(abs(x).^2), 100);
    px = pdf_saddlepoint_approx(y, D, xnd.', varASE, varTherm, sim.verbose);
    if sim.verbose
        for k = 1:length(Pthresh)
            plot(Pthresh(k)*[1 1], [0, 5e4])
        end
        for k = 1:length(Plevels)
            plot(Plevels(k)*[1 1], [0, 10e4])
        end
        for k= 1:length(yd)
            plot(yd(k)*[1 1], [0, 10e4], '--k')
        end
    end
    
    pe2 = 0;
    for k = 1:Nsymb
        if dat(k) == mpam.M-1
            pe2 = pe2 + trapz(y(y < Pthresh(end)), px(y < Pthresh(end), k));
        elseif dat(k) == 0
            pe2 = pe2 + trapz(y(y > Pthresh(1)), px(y > Pthresh(1), k));
        else 
            pe2 = pe2 + trapz(y(y > Pthresh(dat(k) + 1)), px(y > Pthresh(dat(k) + 1), k));
            pe2 = pe2 + trapz(y(y < Pthresh(dat(k))), px(y < Pthresh(dat(k)), k));
        end
    end
    pe2 = real(pe2)/Nsymb;
    
    berpdf = pe2/log2(mpam.M);   
end



    