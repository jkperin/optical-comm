function ber = ber_apd_montecarlo(mpam, tx, fiber, apd, rx, sim)
%% Pre calculations
% Channel response
% Hch does not include transmitter or receiver filter
if isfield(tx, 'modulator')
    Hch = tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay)...
    .*fiber.H(sim.f, tx).*apd.H(sim.f);
%     Hch = tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay)...
%     .*fiber.Hlarge_signal(sim.f, tx).*apd.H(sim.f);
else
    Hch = fiber.H(sim.f, tx).*apd.H(sim.f);
%     Hch = fiber.Hlarge_signal(sim.f, tx).*apd.H(sim.f);
end

link_gain = apd.Gain*apd.R*fiber.link_attenuation(tx.lamb); % Overall link gain

% Ajust levels to desired transmitted power and extinction ratio
mpam = mpam.adjust_levels(tx.Ptx, tx.rexdB);
Pmax = mpam.a(end); % used in the automatic gain control stage

%% Modulated PAM signal
dataTX = randi([0 mpam.M-1], 1, sim.Nsymb); % Random sequence
xt = mpam.mod(dataTX, sim.Mct);
xt(1:sim.Mct*sim.Ndiscard) = 0; % zero sim.Ndiscard first symbols
xt(end-sim.Mct*sim.Ndiscard+1:end) = 0; % zero sim.Ndiscard last symbbols

%% Generate optical signal
[Et, ~] = optical_modulator(xt, tx, sim);

%% Fiber propagation
Et = fiber.linear_propagation(Et, sim.f, tx.lamb);

%% Detect and add noises
yt = apd.detect(Et, sim.fs, 'gaussian', rx.N0);

%% Whitening filter
if sim.WhiteningFilter
    [Hw, yt] = apd.Hwhitening(sim.f, tx.Ptx, rx.N0, yt);
else
    Hw = 1;
end

%% Automatic gain control
% Normalize signal so that highest level is equal to 1
yt = yt/(Pmax*link_gain);
mpam = mpam.norm_levels;

%% Equalization
rx.eq.TrainSeq = dataTX;
[yd, rx.eq] = equalize(rx.eq, yt, Hw.*Hch, mpam, rx, sim);

% Symbols to be discarded in BER calculation
Ndiscard = sim.Ndiscard*[1 1];
if strfind(lower(rx.eq.type), 'adaptive') % must increase discareded symbols to account for filter length, adaptation period, etc
    Ndiscard = Ndiscard + rx.eq.Ndiscard; 
end
ndiscard = [1:Ndiscard(1) sim.Nsymb-Ndiscard(2):sim.Nsymb];
yd(ndiscard) = []; 
dataTX(ndiscard) = [];

%% Demodulate
dataRX = mpam.demod(yd);

%% True BER
[~, ber] = biterr(dataRX, dataTX);

if isfield(sim, 'plots') && sim.plots('Empirical noise pdf')
    % Empirical pdf for a level
    figure(100)
    [nn, xx] = hist(yd(dataTX == 2), 50);
    nn = nn/trapz(xx, nn);
    bar(xx, nn)
    title('Empirical pdf for PAM level 2')
end    

if isfield(sim, 'plots') && sim.plots('Frequency Response')   
    figure(10), hold on, box on    
    leg = {};
    if isfield(tx, 'modulator')
        Hmod = tx.modulator.H(sim.f);
        plot(sim.f/1e9, abs(Hmod).^2)
        leg = [leg 'Modulator'];
    end
        
    if fiber.D(tx.lamb) ~= 0
        Hfib = fiber.H(sim.f, tx);
        plot(sim.f/1e9, abs(Hfib).^2)
        leg = [leg 'Fiber'];
    end
    
    if ~isinf(apd.BW)
        Hapd = apd.H(sim.f);
        plot(sim.f/1e9, abs(Hapd).^2)
        leg = [leg 'APD'];
    end
    
    if sim.WhiteningFilter
        plot(sim.f/1e9, abs(Hw).^2)
        leg = [leg 'Noise Whitening'];
    end
    
    if ~strcmpi(rx.eq.type, 'none')
        plot(sim.f/1e9, abs(rx.eq.Hrx).^2)
        plot(sim.f/1e9, abs(rx.eq.Hff(sim.f/(rx.eq.ros*mpam.Rs))).^2)
        leg = [leg 'Receiver filter', 'Equalizer'];       
    end

    xlabel('Frequency (GHz')
    ylabel('Frequency Response')
    a = axis;
    axis([0 mpam.Rs/1e9 a(3) a(4)])
    legend(leg)
end