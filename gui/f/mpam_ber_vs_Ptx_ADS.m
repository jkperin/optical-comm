function ber_ads = mpam_ber_vs_Ptx_ADS(mpam, tx, fiber1, rx, sim)
addpath ../f % general functions
addpath ../apd
addpath ../apd/f
addpath data/

dBm2Watt = @(x) 1e-3*10.^(x/10);

load SampleData4PAMlong

%% Matlab Simulation
% Simulation parameters
sim.Nsymb = Nsymb; % Number of symbols in montecarlo simulation
sim.Mct = Mct;     % Oversampling ratio to simulate continuous time (must be odd) 
sim.N = N; % number points in 'continuous-time' simulation

% Time and frequency
sim.fs = mpam.Rs*Mct;  % sampling frequency in 'continuous-time'

dt = Ts;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df+0.1).';

sim.t = t;
sim.f = f;

% Transmitter
tx.rexdB = -5;          % exctinction ratio in dB. Defined as Pmin/Pmax
tx_ads = tx;
tx_ads = rmfield(tx_ads, 'modulator');

% Receiver
pin = apd(0, 0, Inf, rx.R, rx.Id);

%% Receiver
if strcmp(rx.filter.type, 'matched')
    % Electric Lowpass Filter
    rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
else
    rx.elefilt = design_filter(rx.filter.type, rx.filterN, rx.filterBw/(sim.fs/2));
end

% Equalization
eq = rx.eq;

% Overall link gain
link_gain = fiber1.link_attenuation(tx.lamb)*pin.R;

% Zeroed samples 
ndiscard = [1:(sim.Ndiscard+eq.Ntrain) (sim.Nsymb-sim.Ndiscard-eq.Ntrain):sim.Nsymb];

%
dataTXref = mpam.demod(xtref(sim.Mct/2+1:sim.Mct:length(xtref)));

% Equalization
eq.TrainSeq = dataTXref;

Ptx = dBm2Watt(tx.PtxdBm);
for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);

	% Ajust levels to desired transmitted power and extinction ratio
    Ptads = Pout.'/mean(Pout)*tx.Ptx;
    
    Ptads([sim.Mct*sim.Ndiscard end-sim.Mct*sim.Ndiscard:end]) = 0; % zero sim.Ndiscard last symbbols

	% Generate optical signal
    [Etads, ~] = optical_modulator(Ptads, tx_ads, sim); % without frequency response, which is already included in ADS simulation

	% Fiber propagation
    [~, Ptads] = fiber1.linear_propagation(Etads, sim.f, tx.lamb);

    % Detect and add noises
    ytads = pin.detect(Ptads, sim.fs, 'gaussian', rx.N0);
  
    % Automatic gain control
    Pmax = 2*tx.Ptx/(1 + 10^(-abs(tx.rexdB)/10));
    ytads = ytads/(Pmax*link_gain);
    mpam.norm_levels;
    
    % Equalization
    ydads = equalize(eq.type, ytads, mpam, tx, fiber1, rx, eq, sim);

    ydads(ndiscard) = [];
    dataTX = dataTXref;
    dataTX(ndiscard) = [];
    
    % Detect ADS signal
    if mpam.M == 2
        ydads = ydads - mean(ytads);
        dataRXads = zeros(size(ydads)).';
        dataRXads(ydads > 0) = 1;
    else        
%         mpam.b = Pthresh.'/mean(Pout)*tx.Ptx;
%         dataRXads = mpam.demod(ydads);
%         mpam.b = btemp;
        dataRXads = mpam.demod(ydads);
    end
    
    % BER
    [~, ber_ads(k)] = biterr(dataRXads, dataTX);
end

% [H, w] = freqz(Heq(end:-1:1), 1);
% figure, hold on
% plot(2*mpam.Rs/1e9*w/(2*pi), abs(H).^2)
% 
% figure, hold on, grid on, box on
% plot(tx.PtxdBm, log10(ber.sim), '-o')
% plot(tx.PtxdBm, log10(ber.ads), '-o')
% plot(tx.PtxdBm, log10(ber.awgn), '-')
% xlabel('Transmitted Power (dBm)')
% ylabel('log_{10}(BER)')
% legend('Matlab Simulation', 'ADS-DFB model', 'AWGN')
% axis([tx.PtxdBm([1 end]) -5 0])



