%% Generate and detect DC-OFDM 
% input:
% ofdm (OFDM parameters)
% tx (transmitter parameters)
% fiber (fiber parameters)
% rx (receiver parameters)
% sim (simulation parameters)

% output:
% ber : counted (including 95%-confidence intervals), and estimated
% Ptx : measured transmitted power

% Function assumes variable number of input arguments just to be compatible
% with previous codes
function [ber, Ptx, Ptx_est] = dc_ofdm(ofdm, tx, fiber, rx, sim)       
       
%% Calculate cyclic prefix ##!! fiber has to be included here
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
Hmod = tx.modulator.H(fc);
Hmod = Hmod.*exp(1j*2*pi*fc*tx.modulator.grpdelay);

% The group delay is removed from the frequency response of the ADC and DAC
Gdac = tx.filter.H(fc/sim.fs);                    
Gadc = rx.filter.H(fc/sim.fs); 

Hfiber = fiber.H(fc, tx);

% Frequency response of the channel at the subcarriers
Gch = K*Gdac.*Hmod.*Hfiber.*rx.R.*Gadc;            

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

%% Adjust signal to get to desired power
rex = 10^(-abs(tx.rexdB)/10);  % defined as Pmin/Pmax    

xmean  = mean(xtc);   % measured

%% Adjust transmitted power
if isfield(tx, 'Ptx') % if Ptx was provided, then scale signal to desired Ptx      
    % Scale and add additional dc bias due to finite extinction ratio
    xtc = xtc*tx.Ptx*(1 - 2*rex)/xmean + 2*tx.Ptx*rex; 
else % just add additional dc bias due to finite extinction ratio
    xtc = xtc + 2*xmean*rex; 
end

%% Apply frequency response of the laser
Et = optical_modulator(xtc, tx, sim);

Ptx = mean(abs(Et).^2);

%% Fiber propagation
Et = fiber.linear_propagation(Et, sim.f, tx.lamb);

%% Direct detection
It = rx.R*abs(Et).^2;

%% Add shot noise
if isfield(sim, 'shot') && sim.shot
    q = 1.60217657e-19;      % electron charge (C)
    Id = 0;                  % dark current

    % Instataneous received power considering only attenuation from the fiber   
    Sshot = 2*q*(It + Id);     % one-sided shot noise PSD

    % Frequency is divided by two because PSD is one-sided
    ws = sqrt(Sshot*sim.fs/2).*randn(size(It));

    It = It + ws;
    
    varshot = mean(Sshot)*sim.fs/2*abs(Gadc).^2/sim.Mct; % referred to after the ADC
else
    varshot = 0;
end

%% Add thermal noise
wt = sqrt(rx.Sth*sim.fs)*randn(sim.N, 1);
It = It + wt;

%% OFDM detection
% If power was readjusted (i.e., Ptx was defined) then calculate scaling
% factor
rex = 10^(tx.rexdB/10);  % defined as Pmin/Pmax     

if isfield(tx, 'Ptx') % if Ptx was provided, then scale signal to desired Ptx
    A = (tx.Ptx*(1 - 2*rex))/(clip_level); % Scale between desired power and actual power
    
    Ptx_est = tx.Ptx;     % Calculate equivalent transmitted power
else
    A = 1;
    
    Ptx_est = clip_level*(1 + 2*rex);
end

ber = ofdm.detect(It, A*Gch, rx, sim);

if isfield(sim, 'RIN') && sim.RIN
    Srin = 10^(tx.RIN/10)*Ptx.^2;
    varrin = Srin*sim.fs*abs(Hfiber.*Gadc).^2/sim.Mct; % referred to after the ADC
else
    varrin = 0;
end
% estimate_ber(A, varthermal, varshot, varrin, Gch)
ber.est = ofdm.estimate_ber(A, rx.Sth*sim.fs*abs(Gadc).^2/sim.Mct, varshot, varrin, Gch);


%% Frequency response
if isfield(sim, 'verbose') && sim.verbose
    figure
    subplot(211), box on, grid on, hold on
    plot(f/1e9, abs(tx.filter.H(f/sim.fs)).^2)
    plot(f/1e9, abs(tx.modulator.H(f)).^2)
    plot(f/1e9, abs(fiber.H(f, tx)).^2)
    plot(f/1e9, abs(rx.filter.H(f/sim.fs)).^2)
    legend('DAC', 'Modulator', 'Fiber', 'Antialiasing (ADC)')
    xlabel('Frequency (GHz)')
    ylabel('|H(f)|^2')

    subplot(212), box on, grid on, hold on
    plot(f/1e9, 20*log10(abs(tx.filter.H(f/sim.fs))))
    plot(f/1e9, 20*log10(abs(tx.modulator.H(f))))
    plot(f/1e9, 20*log10(fiber.H(f, tx)))
    plot(f/1e9, 20*log10(abs(rx.filter.H(f/sim.fs))))
    legend('DAC', 'Modulator', 'Fiber', 'Antialiasing (ADC)')
    xlabel('Frequency (GHz)')
    ylabel('10log_{10}(|H(f)|^2)')
    axis([0 sim.fs/2e9 -40 5])
end




