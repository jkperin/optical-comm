%% Calculate BER for unamplified IM-DD link with APD
% bertail = BER using saddlepoint approximation for the tail probability. 
% Although it's much faster than using the pdf, it can be innacurate due to
% the singularity at the origin.

% bergauss = BER using gaussian approximation

% berpdf (optional) = calculate BER using saddlepoint approximation for the
% pdf, and then calculate tail probability using numerical integration
function [bertail, bergauss, berpdf] = ber_apd_doubly_stochastic(mpam, tx, apd, rx, sim)

Nsymb = mpam.M^sim.L; % number of data symbols 
Nzero = sim.L; % number of zero symbols pad to begin and end of sequnece.
N = sim.Mct*(Nsymb + 2*Nzero); % total number of points 

% Frequency
df = 1/N;
f = (-0.5:df:0.5-df).';

% Generate optical signal
try 
    dataTX = debruijn_generator(mpam.M, sim.L).' - 1;
catch e % can fail for some states of the random number generator
    warning(e.message)
    load rng_state  % state that works
    rng(rng_state);
    dataTX = debruijn_generator(mpam.M, sim.L).' - 1;
    rng('shuffle')
end

% Generate symbols
Pd = mpam.a(gray2bin(dataTX, 'pam', mpam.M) + 1);

% modulate
Pt = reshape(kron(Pd, mpam.pshape(1:sim.Mct)).', [], 1);

% Rescale to desired power
if strcmp(mpam.level_spacing, 'uniform') % adjust levels to include extinction ratio penalty
    Pmin = mean(Pt)/(10^(abs(tx.rex)/10)-1); % minimum power 
    Plevels = (mpam.a + Pmin)*apd.R*apd.Gain*tx.Ptx/(mean(Pt) + Pmin); % levels at the receiver
    Pthresh = (mpam.b + Pmin)*apd.R*apd.Gain*tx.Ptx/(mean(Pt) + Pmin); % decision thresholds at the receiver
    Pt = (Pt + Pmin)*apd.R*apd.Gain*tx.Ptx/(mean(Pt) + Pmin); % after rescaling E(Pt) = tx.Ptx
elseif strcmp(mpam.level_spacing, 'nonuniform') % already includes extinction ratio penalty, so just scale
    Plevels = mpam.a*apd.R*apd.Gain*tx.Ptx/mean(Pt);
    Pthresh = mpam.b*apd.R*apd.Gain*tx.Ptx/mean(Pt);
    Pt = Pt*apd.R*apd.Gain*tx.Ptx/mean(Pt);
end  

yt = [zeros(Nzero*sim.Mct, 1); Pt; zeros(Nzero*sim.Mct, 1)]; % zero pad

% Electric low-pass filter
yt = real(ifft(fft(yt).*ifftshift(rx.elefilt.H(f))));

% Discard zeros and downsample
ix = sim.Mct*Nzero+sim.Mct/2:sim.Mct:N-sim.Mct*Nzero;
yd = yt(ix);

if sim.verbose
    figure(102), hold on
    plot([zeros(Nzero*sim.Mct, 1); Pt; zeros(Nzero*sim.Mct, 1)])
    plot(ix, yd, 'o')
    legend('Transmitted power', 'KL-SE Fourier')
end

%%
varTherm = rx.N0*rx.elefilt.noisebw(sim.fs)/2; % variance of thermal noise

pe = 0;
pe_gauss = 0;
dat = gray2bin(dataTX, 'pam', mpam.M);
for k = 1:Nsymb
    mu = yd(k);
    varShot = 2*apd.q*apd.Gain^2*apd.Fa*(apd.R*yd(k)/apd.Gain + apd.Id)*rx.elefilt.noisebw(sim.fs)/2; % Agrawal 4.4.17 (4th edition)
    sig = sqrt(varTherm + varShot);
  
    lambda = (apd.R*yd(k) + apd.Id)*(sim.Mct/sim.fs)/apd.q; 
    
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