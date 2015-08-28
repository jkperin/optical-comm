%% Compare performance of M-PAM in 3 different scenarios
% 1. PIN receiver, no amplifier
% 2. SOA & PIN receiver
% 3. APD
clear, clc, close all

addpath ../f % general functions
addpath ../soa
addpath ../soa/f
addpath ../apd
addpath ../apd/f

%% Simulation parameters
sim.Nsymb = 2^16; % Number of symbols in montecarlo simulation
sim.Mct = 15;     % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
sim.Me = 16; % number of used eigenvalues
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

sim.polarizer = false;
sim.shot = true; % include shot noise in montecarlo simulation (always included for pin and apd case)
sim.RIN = true; % include RIN noise in montecarlo simulation
sim.verbose = false; % show stuff

%% M-PAM
mpam = PAM(4, 100e9, 'optimized', @(n) double(n >= 0 & n < sim.Mct));

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

%% Transmitter
switch mpam.M
    case 4
        tx.PtxdBm = -24:1:-8;
    case 8
        tx.PtxdBm = -22:2:-4;
    case 16
       tx.PtxdBm = -18:2:-2;
end
   
tx.lamb = 1310e-9; % wavelength
tx.alpha = 0; % chirp parameter
tx.RIN = -150;  % dB/Hz
tx.rexdB = -15;  % extinction ratio in dB. Defined as Pmin/Pmax

%% Fiber
fiber = fiber(); % fiber(L, att(lamb), D(lamb))

%% Receiver
rx.N0 = (30e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
% rx.matchedfilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
% rx.elefilt = design_filter('matched', @(t) conv(mpam.pshape(t), 1/sim.fs*tx.modulator.h(t/sim.fs), 'full') , 1/sim.Mct);
% Optical Bandpass Filter
rx.optfilt = design_filter('fbg', 0, 200e9/(sim.fs/2));

%% Equalization
rx.eq.type = 'Fixed TD-SR-LE';
% rx.eq.ros = 2;
rx.eq.Ntaps = 31;
% rx.eq.Ntrain = 2e3;
% rx.eq.mu = 1e-2;

% KLSE Fourier Series Expansion (done here because depends only on filters
% frequency response)
% klse_fourier(rx, sim, N, Hdisp)
[rx.U_fourier, rx.D_fourier, rx.Fmax_fourier] = klse_fourier(rx, sim, sim.Mct*(mpam.M^sim.L + 2*sim.L));

%% PIN
% (GaindB, ka, GainBW, R, Id) 
pin = apd(0, 0, Inf, rx.R, rx.Id);

%% APD 
% (GaindB, ka, GainBW, R, Id) 
% Finite Gain x BW
apd_fix = apd(10*log10(5), 0.09, 340e9, rx.R, rx.Id); 

% % Infinite Gain x BW
% apd_opt = apd(7, 0.09, Inf, rx.R, rx.Id);
% % % Optimized Gain
% sim.awgn = true;
% apd_opt.optimize_gain(mpam, tx, fiber, rx, sim);
% sim = rmfield(sim, 'awgn');

%% SOA
% soa(GaindB, NF, lambda, maxGaindB)
soa = soa(20, 7, 1310e-9, 20); 

%
Fc = (20:5:50)*1e9;

for k = 1:length(Fc)
%     % Modulator frequency response
    tx.modulator.fc = Fc(k); % modulator cut off frequency
    tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
    tx.modulator.h = @(t) [0*t(t < 0) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0))];
    tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds
   
    [~, eq] = equalize(rx.eq.type, [], mpam, tx, fiber, rx, sim);
          
    % BER
    disp('BER with SOA')
    ber.soa(k) = soa_ber(mpam, tx, fiber, soa, rx, sim);
    disp('BER with APD with fix gain')
    ber.apd_fix(k) = apd_ber(mpam, tx, fiber, apd_fix, rx, sim);
%     disp('BER with APD with optimized gain')
%     ber.apd_opt(k) = apd_ber(mpam, tx, fiber, apd_opt, rx, sim);
    disp('BER with PIN')
    ber.pin(k) = apd_ber(mpam, tx, fiber, pin, rx, sim);
    
    % Plot
    figure, hold on, grid on, box on
    plot(tx.PtxdBm, log10(ber.soa(k).est), '-b')
    plot(tx.PtxdBm, log10(ber.apd_fix(k).gauss), '-r')
%     plot(tx.PtxdBm, log10(ber.apd_opt(k).gauss), '-m')
    plot(tx.PtxdBm, log10(ber.pin(k).gauss), '-k')

    plot(tx.PtxdBm, log10(ber.soa(k).count), '--ob')
    plot(tx.PtxdBm, log10(ber.apd_fix(k).count), '--or')
%     plot(tx.PtxdBm, log10(ber.apd_opt(k).count), '--om')
    plot(tx.PtxdBm, log10(ber.pin(k).count), '--ok')

%     plot(tx.PtxdBm, log10(ber.soa(k).gauss), '--b')

   plot(tx.PtxdBm, log10(ber.soa(k).awgn_ne), ':b')
    plot(tx.PtxdBm, log10(ber.apd_fix(k).awgn_ne), ':r')
%     plot(tx.PtxdBm, log10(ber.apd_opt(k).awgn_ne), ':m')
    plot(tx.PtxdBm, log10(ber.pin(k).awgn_ne), ':k')
       
    Preq.pin(k) = interp1(log10(ber.pin(k).gauss), tx.PtxdBm, log10(sim.BERtarget), 'spline');
    Preq.apd_fix(k) = interp1(log10(ber.apd_fix(k).gauss), tx.PtxdBm, log10(sim.BERtarget), 'spline');
%     Preq.apd_opt(k) = interp1(log10(ber.apd_opt(k).gauss), tx.PtxdBm, log10(sim.BERtarget), 'spline');
    Preq.soa(k) = interp1(log10(ber.soa(k).awgn_ne), tx.PtxdBm, log10(sim.BERtarget), 'spline');
    
%     % Get only valid points
%     log_ber_soa = log10(ber.soa(k).count);
%     ix_soa = ~(isinf(log_ber_soa) | isnan(log_ber_soa));
%     PdBm = tx.PtxdBm(ix_soa);
%     
%     Pfit = interp1(log_ber_soa(ix_soa), PdBm, -6:0) ;
%     plot(Pfit, -6:0, 'g')
%     
%     Preq.soa(k) = interp1(log_ber_soa(ix_soa), PdBm, log10(sim.BERtarget));

    plot(Preq.pin(k), log10(sim.BERtarget), 'dg')
    plot(Preq.apd_fix(k), log10(sim.BERtarget), 'dg')
%     plot(Preq.apd_opt(k), log10(sim.BERtarget), 'dg')
    plot(Preq.soa(k), log10(sim.BERtarget), 'dg')
    
    xlabel('Received Power (dBm)')
    ylabel('log(BER)')
    legend('SOA', 'APD fixed gain', 'APD optimized gain', 'PIN', 'Location', 'SouthWest')
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -8 0])
    set(gca, 'xtick', tx.PtxdBm)
    title(sprintf('Modulator BW = %d GHz', Fc(k)/1e9))
    
end

figure, hold on, box on, grid on
% plot(Fc/1e9, Preq.pin-Preq.apd_opt, 'b')
plot(Fc/1e9, Preq.pin-Preq.apd_fix, 'g')
plot(Fc/1e9, Preq.pin-Preq.soa, 'r')
xlabel('Modulator Bandwidth (GHz)')
ylabel('Margin Improvement (dB)')
legend('APD fix', 'SOA')

Preq.pin_eqs = [-9.857216456013484,-10.609901105249634,-11.066972526647389,-11.362463326153650,-11.565995005086206,-11.712651990847231,-11.822208550113094];
Preq.apd_fix_eqs = [-15.878872818617255,-16.805442413152400,-17.351791802193482,-17.701175021238033,-17.938687086993948,-18.107960179223504,-18.233696828430034];
% Preq.soa_eps = [-15.825386703118062,-17.312834244944270,-18.213914547752620,-18.797955189186148,-19.198430060160150,-19.486094187948530,-19.700898375907624];
Preq.soa_eqs = [-16.5, -18.75, -19, -19.4, -19.5, -20.8, -20.8];
% SOA has to be measured meanually with ber.count

Preq.pin_opt = [-9.874932457075868,-10.625855517152765,-11.082025617775482,-11.377068508699681,-11.580352979174755,-11.726827874342657,-11.836234506673815];
Preq.apd_fix_opt = [-16.210564591028240,-17.088290849318230,-17.607315289928770,-17.940325684073503,-18.166274980839532,-18.327747609965960,-18.448066504987230];
% Preq.soa_opt = [-15.895049803922603,-17.336015557482970,-18.209518282923348,-18.775285555490242,-19.162921007419257,-19.441106735827898,-19.648673749931230];
Preq.soa_opt = [-13.5, -16.5, -18.5, -19.8, -20, -21, -21.5];
% % 
% % 
figure, hold on, box on, grid on
plot(Fc/1e9, Preq.pin_eqs-Preq.apd_fix_eqs, 'b')
% plot(Fc/1e9, Preq.pin_eqs-Preq.apd_opt_eqs, 'g')
plot(Fc/1e9, Preq.pin_eqs-Preq.soa_eqs, 'r')

plot(Fc/1e9, Preq.pin_eqs-Preq.apd_fix_opt, '--b')
% plot(Fc/1e9, Preq.pin_eqs-Preq.apd_opt_opt, '--g')
plot(Fc/1e9, Preq.pin_eqs-Preq.soa_opt, '--r')
xlabel('Modulator Bandwidth (GHz)')
ylabel('Margin Improvement (dB)')
% axis([20 50 0 8])
legend('APD Gain = 7 dB', 'SOA', 'Location', 'SouthEast')


 %% 
% Preq.pin_eqs = [-9.856962921183065,-10.609844934281513,-11.066957417744230,-11.362458826005494,-11.565993507847503,-11.712651912409520,-11.822208260844207];
% Preq.apd_fix_eqs = [-15.883321186494720,-16.812194481362926,-17.359278932653800,-17.708982078139800,-17.946641021696745,-18.115997091258144,-18.241802907673573];
% Preq.apd_opt_eqs = [-16.630628660627018,-17.653792716607150,-18.252406291960300,-18.632768155727180,-18.890531527110806,-19.074010866580736,-19.210084312694537];
% Preq.soa_eqs = [-14.829286758052342,-16.317549588264950,-17.219277561398325,-17.803814708619860,-18.204657129491090,-18.492609163326538,-18.707646822333988];
% 
% Preq.pin_opt = [-9.874910852208039,-10.625864595588867,-11.082029722036864,-11.377070026835572,-11.580353496330398,-11.726828647403268,-11.836234525239737];
% Preq.apd_fix_opt = [-16.217949810873815,-17.096719315151980,-17.616079885594190,-17.949212431266640,-18.175216527350810,-18.336753212562430,-18.457122569406152];
% Preq.apd_opt_opt = [-17.205437853681993,-18.149278585784618,-18.701710938765570,-19.053496947743973,-19.291414164267508,-19.461236990110184,-19.587500808245515];
% % Preq.soa_opt = [-15.655864144465879,-16.690378567965380,-17.413535263709100,-17.910509730932580,-18.260898134784327,-18.516511946023485,-18.709143859942827];
% % Preq.soa_opt = [-14,-16,-17,-19.5,-19,-20,-20.5];
% Preq.soa_opt = [-16.6093  -18.0975  -18.9993  -19.5838  -19.9847  -20.2726  -20.4876];
% %opt (50, -20.5), (45, -20), (