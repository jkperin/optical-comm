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

% Channel Response
Ptx = design_filter('matched', mpam.pshape, 1/sim.Mct); % transmitted pulse shape

% Hch does not include receiver filter
if isfield(tx, 'modulator')
    Hch = Ptx.H(sim.f/sim.fs).*tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay)...
    .*fiber.H(sim.f, tx).*apd.H(sim.f);
else
    Hch = Ptx.H(sim.f/sim.fs).*fiber.H(sim.f, tx).*apd.H(sim.f);
end

link_gain = apd.Gain*apd.R*fiber.link_attenuation(tx.lamb); % Overall link gain. Equivalent to Hch(0)

% Ajust levels to desired transmitted power and extinction ratio
mpam = mpam.adjust_levels(tx.Ptx, tx.rexdB);
Pmax = mpam.a(end); % used in the automatic gain control stage

% Modulated PAM signal
dataTX = debruijn_sequence(mpam.M, sim.L).'; % de Bruijin sequence
xt = mpam.mod(dataTX, sim.Mct);
xt = [zeros(sim.Mct*Nzero, 1); xt; zeros(sim.Mct*Nzero, 1)]; % zero pad

% Generate optical signal
if ~isfield(sim, 'RIN')
    sim.RIN = false;
end

RIN = sim.RIN;
sim.RIN = false; % RIN is not modeled here since number of samples is not high enough to get accurate statistics
[Et, ~] = optical_modulator(xt, tx, sim);

% Fiber propagation
[~, Pt] = fiber.linear_propagation(Et, sim.f, tx.lamb);

% Direct detect
yt = apd.detect(Pt, sim.fs, 'no noise');

% Automatic gain control
% Pmax = 2*tx.Ptx/(1 + 10^(-abs(tx.rexdB)/10)); % calculated from mpam.a
yt = yt/(Pmax*link_gain); % just refer power values back to transmitter
mpam = mpam.norm_levels;

%% Equalization
if isfield(rx, 'eq') && (isfield(tx, 'modulator') || ~isinf(apd.BW))
    rx.eq.type = strrep(rx.eq.type, 'Adaptive', 'Fixed'); % replace adaptive for fixed
    % Note: in this simulation there aren't many symbols to determine the
    % adaptive equalizer
else % otherwise only filter using rx.elefilt
    rx.eq.type = 'None';
end

% Equalizer
[yd, rx.eq] = equalize(rx.eq, yt, Hch, mpam, rx, sim);

% Symbols to be discard in BER calculation
yd = yd(Nzero+1:end-Nzero);

%% Detection
Pthresh = mpam.b; % refer decision thresholds to receiver

% Noise bandwidth
Heq = freqz(rx.eq.num, rx.eq.den, sim.f, mpam.Rs);
Htot = Hch.*Heq; 
Df  = (1/(abs(Htot(sim.f==0)).^2))*trapz(sim.f, abs(Htot).^2)/2; % matched filter noise BW
Dfshot = (1/(abs(apd.H(0).*Htot(sim.f==0)).^2))*trapz(sim.f, abs(apd.H(sim.f).*Htot).^2)/2;
DfRIN = (1/(abs(fiber.H(0, tx).*apd.H(0).*Htot(sim.f==0)).^2))*trapz(sim.f, abs(fiber.H(sim.f, tx).*apd.H(sim.f).*Htot).^2)/2;

% Variance of thermal noise
varTherm = rx.N0*Df; 

if RIN
    varRIN =  @(Pnorm) 10^(tx.RIN/10)*(Pmax*link_gain*Pnorm).^2*DfRIN;
else
    varRIN = @(Pnorm) 0;
end

pe = 0;
pe_gauss = 0;
dat = gray2bin(dataTX, 'pam', mpam.M);
for k = 1:Nsymb
    mu = yd(k);
    varShot = apd.varShot(Pmax*link_gain*yd(k)/apd.Gain, Dfshot); 
    % Note: apd.BW*pi/2 is the noise BW of the APD frequency response
    sig = 1/(Pmax*link_gain)*sqrt(varTherm + varShot + varRIN(yd(k)));
     
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