%% Generate and detect DC-OFDM 
% input:
% ofdm (OFDM parameters)
% tx (transmitter parameters)
% fiber (fiber parameters)
% rx (receiver parameters)
% sim (simulation parameters)

% output:
% ber: counted (including 95%-confidence intervals), and estimated

% Function assumes variable number of input arguments just to be compatible
% with previous codes
function ber = dc_ofdm(ofdm, tx, fiber, rx, sim)       
       
%% Inserts cyclic prefix ##!! fiber has to be included here
ofdm.cyclic_prefix(tx, rx, sim);

%% Time and frequency scales
sim.fs = sim.Mct*ofdm.fs;                                     % sampling frequency to emulate continuous time (Hz)  

sim.N = sim.Mct*(ofdm.Nc + ofdm.Npre_os)*sim.Nsymb;           % total number of points simulated in continuous time
         
dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;  

%% Frequency responses
fc = ofdm.fc; % subcarrier frequencies

K = 1 - 2*qfunc(tx.rclip);  % Amplitude attenuation due to clipping = (1-Q(r))

% Remove group delay of modulator frequency response
Hmod = tx.kappa*tx.modulator.H(fc);
Hmod = Hmod.*exp(1j*2*pi*fc*tx.modulator.grpdelay);

% The group delay is removed from the frequency response of the ADC and DAC
Gdac = tx.filter.H(fc/sim.fs);                    
Gadc = rx.filter.H(fc/sim.fs); 

Hfiber = fiber.Hfiber(fc, tx);

% Frequency response of the channel at the subcarriers
Gch = K*Gdac.*tx.kappa.*Hmod.*Hfiber.*rx.R.*Gadc;            

%% Power allocation
ofdm.power_allocation(tx, fiber, rx, sim);

% Signal std at transmitter and receiver
sigtx = sqrt(2*sum(ofdm.Pn));

%% Define clipping levlel
% (i) full dc bias: sigtx*rclip
% (ii) reduced dc bias: sigtxdac*rclip
if isfield(sim, 'full_dc') && sim.full_dc
    clip_level = tx.rclip*sigtx;
else
    sim.full_dc = false;
    clip_level = tx.rclip*sqrt(sum(2*ofdm.Pn.*abs(Gdac).^2));
end

% Check if power is at reasonable levels
assert(sigtx < 1e-3, 'Power is too high: %.2G (mW)\nCannot improve performance by increasing power.\n', 1e3*sqrt(2*sum(ofdm.Pn)))

%% Generate OFDM signal
xt = ofdm.generate_signal(tx, sim);

%% Clipping
xtc = xt;
xtc(xt < -clip_level) = -clip_level;
xtc(xt > clip_level) = clip_level;

%% Add dc bias (dc_bias determines whether it's added full dc bias or clipped dc bias)
xtc = xtc + clip_level;

%% Apply frequency response of the laser
Et = optical_modulator(xtc, tx, sim);

%% Fiber propagation
Et = fiber.linear_propagation(Et, sim.f, tx.lamb);

%% Direct detection
It = rx.R*abs(Et).^2;

%% Add shot noise
if sim.shot
    q = 1.60217657e-19;      % electron charge (C)
    Id = 0;                  % dark current

    % Instataneous received power considering only attenuation from the fiber   
    Sshot = 2*q*(It + Id);     % one-sided shot noise PSD

    % Frequency is divided by two because PSD is one-sided
    ws = sqrt(Sshot*sim.fs/2).*randn(size(It));

    It = It + ws;
end

%% Add thermal noise
wt = sqrt(rx.Sth*sim.fs)*randn(sim.N, 1);
It = It + wt;

%% OFDM detection
% If power was readjusted (i.e., Ptx was defined) then calculate scaling
% factor
if isfield(tx, 'Ptx') % if Ptx was provided, then scale signal to desired Ptx
    % Calculate equivalent transmitted power
    Ptx  = tx.kappa*clip_level;   % = Ptx if Hmod(f=0) = 1

    rex = 10^(tx.rexdB/10);  % defined as Pmin/Pmax      

    % Scale and add additional dc bias corresponding to finite
    % extinction ratio
    A = (tx.Ptx*(1 - 2*rex))/Ptx;
else
    A = 1;
end

ber = ofdm.detect(It, A*Gch, rx, sim);


%% Frequency response
if sim.verbose
    figure
    subplot(211), box on, grid on, hold on
    plot(f/1e9, abs(tx.filter.H(f/sim.fs)).^2)
    plot(f/1e9, abs(tx.modulator.H(f)).^2)
    plot(f/1e9, abs(fiber.Hfiber(f, tx)).^2)
    plot(f/1e9, abs(rx.filter.H(f/sim.fs)).^2)
    legend('DAC', 'Modulator', 'Fiber', 'Antialiasing (ADC)')
    xlabel('Frequency (GHz)')
    ylabel('|H(f)|^2')

    subplot(212), box on, grid on, hold on
    plot(f/1e9, 20*log10(abs(tx.filter.H(f/sim.fs))))
    plot(f/1e9, 20*log10(abs(tx.modulator.H(f))))
    plot(f/1e9, 20*log10(fiber.Hfiber(f, tx)))
    plot(f/1e9, 20*log10(abs(rx.filter.H(f/sim.fs))))
    legend('DAC', 'Modulator', 'Fiber', 'Antialiasing (ADC)')
    xlabel('Frequency (GHz)')
    ylabel('10log_{10}(|H(f)|^2)')
    axis([0 sim.fs/2e9 -40 5])
end




