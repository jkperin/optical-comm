%% Calculate BER for amplified IM-DD link
%% Shot and RIN are not included even if sim.RIN and sim.shot == true. 
%% Signal-spontaenous noise is assumed to be dominant followed by thermal noise
% bertail = BER using saddlepoint approximation for the tail probability. 
% Although it's much faster than using the pdf, it can be innacurate due to
% the singularity at the origin.

% bergauss = BER using gaussian approximation

% berpdf (optional) = calculate BER using saddlepoint approximation for the
% pdf, and then calculate tail probability using numerical integration
function [bertail, bergauss, berpdf] = ber_soa_klse_fourier(mpam, tx, fiber, soa, rx, sim)

% From klse_fourier.m
U = rx.U_fourier;
D = rx.D_fourier;
Fmax = rx.Fmax_fourier;

if isfield(sim, 'polarizer') && ~sim.polarizer
    Npol = 2;     % number of noise polarizations
else % by default assumes that polarizer is being use so Npol = 1.
    Npol = 1;
end
Nsymb = mpam.M^sim.L; % number of data symbols 
Nzero = sim.L; % number of zero symbols padded to begin and end of sequnece.
N = sim.Mct*(Nsymb + 2*Nzero); % total number of points 

% Frequency
df = 1/N;
f = (-0.5:df:0.5-df).';
sim.f = f*sim.fs; % redefine frequency to be used in optical_modulator.m and fiber propagation

% Overall link gain
link_gain = soa.Gain*fiber.link_attenuation(tx.lamb)*rx.R;

% Ajust levels to desired transmitted power and extinction ratio
mpam.adjust_levels(tx.Ptx, tx.rexdB);

%% Modulated PAM signal
dataTX = debruijn_sequence(mpam.M, sim.L).'; % de Bruijin sequence
xt = mpam.mod(dataTX, sim.Mct);
xt = [zeros(sim.Mct*Nzero, 1); xt; zeros(sim.Mct*Nzero, 1)]; % zero pad

%% Generate optical signal
RIN = sim.RIN;
sim.RIN = false; % RIN is not modeled here since number of samples is not high enough to get accurate statistics
[Et, Pt] = optical_modulator(xt, tx, sim);

% Fiber propagation moved to KLSE Fourier

%% Signal after amplifier
x = sqrt(soa.Gain)*Et;

% Fourier series coefficients for periodic extension of x 
xn = fftshift(fft(x))/length(x); % divided by period since it's Fourier series
xn = xn(abs(f) <= Fmax);

fm = f(abs(f) <= Fmax);

% Calculate output
xk = zeros(length(fm), N);
for k = 1:length(x)
    xk(:, k) = (xn.*exp(1j*2*pi*fm*(k-1)));
end

% Used to calculate non-centrality parameter of Chi-Squared distributions
% ck(i, j) = ith chi-square distribution at jth time sample. |ck(i, j)|^2
% is the non-centrality parameter of that chi-square distribution
ck = U'*xk; 

% Discard zeros and downsample
ix = sim.Mct*Nzero+1+(sim.Mct-1)/2:sim.Mct:N-sim.Mct*Nzero; % sampling points
ck = ck(:, ix);
yd = abs(x(ix)).^2;

%% Detection
Pthresh = mpam.b*link_gain; % refer decision thresholds to receiver

varASE = (soa.N0*sim.fs/2)/N; % variance in each dimension of the circularly symmetric Gaussian noise
% Note: factor 1/N appears because in this derivation we assume signal 
% component is periodic with period N.
% Note: soa.N0 is single-sided baseband equivalent of ASE PSD 

% Noise bandwidth
Df  = rx.elefilt.noisebw(sim.fs)/2;
% Variance of thermal noise
varTherm = rx.N0*Df; 

% Variance of signal dependent noise
if sim.shot
    varShot = @(Plevel) 2*1.60217657e-19*(rx.R*Plevel + rx.Id)*Df;
else
    varShot = @(Plevel) 0;
end

if RIN
    varRIN =  @(Plevel) 10^(tx.RIN/10)*Plevel.^2*Df;
else
    varRIN = @(Plevel) 0;
end

varSigDependent = @(Plevel) varShot(Plevel) + varRIN(Plevel);

%% BER calculation using saddle-point approximation and Gaussian approximation
pe = 0;
pe_gauss = 0;
dat = gray2bin(dataTX, 'pam', mpam.M);
for k = 1:Nsymb
    alphan = D.*abs(ck(:, k)).^2;
    betan = D*varASE;
    mu = sum(alphan + Npol*betan); % mean of output random variable
    sig = sqrt(sum(betan.*(2*alphan + Npol*betan)) + varTherm + varSigDependent(yd(k))); % std of output random variable (including thermal noise)
    
    if dat(k) == mpam.M-1
        pe = pe + tail_saddlepoint_approx(Pthresh(end), D, ck(:, k), varASE, varTherm + varSigDependent(yd(k)), 'left', Npol);
        
        % Error probability using Gaussian approximation
        pe_gauss = pe_gauss + qfunc((mu-Pthresh(end))/sig);
    elseif dat(k) == 0;
        pe = pe + tail_saddlepoint_approx(Pthresh(1), D, ck(:, k), varASE, varTherm + varSigDependent(yd(k)), 'right', Npol);
        
        % Error probability using Gaussian approximation
        pe_gauss = pe_gauss + qfunc((Pthresh(1)-mu)/sig);
    else 
        pe = pe + tail_saddlepoint_approx(Pthresh(dat(k) + 1), D, ck(:, k), varASE, varTherm + varSigDependent(yd(k)), 'right', Npol);
        pe = pe + tail_saddlepoint_approx(Pthresh(dat(k)), D, ck(:, k), varASE, varTherm + varSigDependent(yd(k)), 'left', Npol);
        
        % Error probability using Gaussian approximation
        pe_gauss = pe_gauss + qfunc((Pthresh(dat(k) + 1) - mu)/sig);
        pe_gauss = pe_gauss + qfunc((mu - Pthresh(dat(k)))/sig);
    end
end

pe = real(pe)/Nsymb;
pe_gauss = real(pe_gauss)/Nsymb;

bertail = pe/log2(mpam.M);
bergauss = pe_gauss/log2(mpam.M);

if sim.verbose
    A = U*diag(D)*U';
    yk = zeros(length(x), 1); % output signal (noiseless)
    for k = 1:length(x)
        yk(k) = real(xk(:, k)'*A*xk(:, k));
    end
    
    yd = yk(ix);
    
    figure(102), hold on
    plot(link_gain*Pt)
    plot(yk, '-k')
    plot(ix, yd, 'o')
    legend('Transmitted power', 'KL-SE Fourier')
end

%% Calculate pdfs using saddlepoint approximation
if nargout == 3
    y = linspace(0, 2*max(abs(x).^2), 100);
    px = pdf_saddlepoint_approx(y, D, ck, varASE, varTherm, Npol, sim.verbose);
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
