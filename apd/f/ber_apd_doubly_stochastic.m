%% Calculate BER for unamplified IM-DD link with APD
% bertail = BER using saddlepoint approximation for the tail probability. 
% Although it's much faster than using the pdf, it can be innacurate due to
% the singularity at the origin.

% bergauss = BER using gaussian approximation

% berpdf (optional) = calculate BER using saddlepoint approximation for the
% pdf, and then calculate tail probability using numerical integration
function [bertail, bergauss, berpdf] = ber_apd_doubly_stochastic(mpam, tx, fiber, apd, rx, sim)

Nsymb = mpam.M^sim.L; % number of data symbols 
% number of zero symbols pad to begin and end of sequnece.
if isfield(rx, 'eq') && isfield(rx.eq, 'Ntaps')
    Nzero = max(rx.eq.Ntaps, sim.L);
else
    Nzero = sim.L;
end
N = sim.Mct*(Nsymb + 2*Nzero); % total number of points 

% Frequency
df = 1/N;
f = (-0.5:df:0.5-df).';
sim.f = f*sim.fs; % redefine frequency to be used in optical_modulator.m

% Overall link gain
link_gain = apd.Gain*fiber.link_attenuation(tx.lamb)*apd.R;

% Ajust levels to desired transmitted power and extinction ratio
mpam.adjust_levels(tx.Ptx, tx.rexdB);
Pmax = mpam.a(end); % used in the automatic gain control stage

% Modulated PAM signal
dataTX = debruijn_sequence(mpam.M, sim.L).'; % de Bruijin sequence
xt = mpam.mod(dataTX, sim.Mct);
xt = [zeros(sim.Mct*Nzero, 1); xt; zeros(sim.Mct*Nzero, 1)]; % zero pad

% Generate optical signal
RIN = sim.RIN;
sim.RIN = false; % RIN is not modeled here since number of samples is not high enough to get accurate statistics
[Et, ~] = optical_modulator(xt, tx, sim);

% Fiber propagation
[~, Pt] = fiber.linear_propagation(Et, sim.f, tx.lamb);

% Direct detect
yt = apd.R*apd.Gain*Pt;

% Automatic gain control
% Pmax = 2*tx.Ptx/(1 + 10^(-abs(tx.rexdB)/10)); % calculated from mpam.a
yt = yt/(Pmax*link_gain); % just refer power values back to transmitter
mpam.norm_levels;

%% Equalization
if isfield(rx, 'eq')
    rx.eq.type = strrep(rx.eq.type, 'Adaptive', 'Fixed'); % replace adaptive for fixed
    % Note: in this simulation there aren't many symbols to determine the
    % adaptive equalizer
else % otherwise only filter using rx.elefilt
    rx.eq.type = 'None';
end

% Equalizer
[yd, rx.eq] = equalize(rx.eq.type, yt, mpam, tx, fiber, rx, sim);

% Symbols to be discard in BER calculation
yd = yd(Nzero+1:end-Nzero);

% Ensures signal is non-negative
yd(yd < 0) = 0; % !! why this here?

%% Detection
Pthresh = mpam.b; % refer decision thresholds to receiver

% Noise bandwidth
Df  = rx.elefilt.noisebw(sim.fs)/2;

% Variance of thermal noise
varTherm = rx.N0*Df; 

if RIN
    varRIN =  @(Pnorm) 10^(tx.RIN/10)*(Pmax*link_gain*Pnorm).^2*Df;
else
    varRIN = @(Pnorm) 0;
end

pe = 0;
pe_gauss = 0;
dat = gray2bin(dataTX, 'pam', mpam.M);
for k = 1:Nsymb
    mu = yd(k);
    varShot = apd.varShot(Pmax*link_gain*yd(k)/apd.Gain, Df);
    sig = 1/(Pmax*link_gain)*sqrt(rx.eq.Kne)*sqrt(varTherm + varShot + varRIN(yd(k)));
     
    if dat(k) == mpam.M-1
%         pe = pe + apd.tail_saddlepoint_approx(Pthresh(end), lambda, sim.fs/sim.Mct, rx.N0, 'left');
        
        % Error probability using Gaussian approximation
        pe_gauss = pe_gauss + qfunc((mu-Pthresh(end))/sig);
    elseif dat(k) == 0;
%         pe = pe + apd.tail_saddlepoint_approx(Pthresh(1), lambda, sim.fs/sim.Mct, rx.N0, 'right');
        
        % Error probability using Gaussian approximation
        pe_gauss = pe_gauss + qfunc((Pthresh(1)-mu)/sig);
    else 
%         pe = pe + apd.tail_saddlepoint_approx(Pthresh(dat(k) + 1), lambda, sim.fs/sim.Mct, rx.N0, 'right');
%         pe = pe + apd.tail_saddlepoint_approx(Pthresh(dat(k)), lambda, sim.fs/sim.Mct, rx.N0, 'left');
        
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
    figure(102), hold on
    plot(link_gain*xt)
    plot(ix, yd, 'o')
    legend('Transmitted power', 'Sampled received signal')
end

%% Calculate pdfs using saddlepoint approximation
if nargout == 3
    berpdf = 0;
%     y = linspace(0, 2*max(abs(x).^2), 100);
%     px = pdf_saddlepoint_approx(y, D, ck, varASE, varTherm, sim.verbose);
%     if sim.verbose
%         for k = 1:length(Pthresh)
%             plot(Pthresh(k)*[1 1], [0, 5e4])
%         end
%         for k = 1:length(Plevels)
%             plot(Plevels(k)*[1 1], [0, 10e4])
%         end
%         for k= 1:length(yd)
%             plot(yd(k)*[1 1], [0, 10e4], '--k')
%         end
%     end
%     
%     pe2 = 0;
%     for k = 1:Nsymb
%         if dat(k) == mpam.M-1
%             pe2 = pe2 + trapz(y(y < Pthresh(end)), px(y < Pthresh(end), k));
%         elseif dat(k) == 0
%             pe2 = pe2 + trapz(y(y > Pthresh(1)), px(y > Pthresh(1), k));
%         else 
%             pe2 = pe2 + trapz(y(y > Pthresh(dat(k) + 1)), px(y > Pthresh(dat(k) + 1), k));
%             pe2 = pe2 + trapz(y(y < Pthresh(dat(k))), px(y < Pthresh(dat(k)), k));
%         end
%     end
%     pe2 = real(pe2)/Nsymb;
%     
%     berpdf = pe2/log2(mpam.M);   
end