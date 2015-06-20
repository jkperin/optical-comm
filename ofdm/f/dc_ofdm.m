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

%% Define and redefine some variables
Nc = ofdm.Nc; 

% The group delay is removed from the frequency response of the ADC and DAC
Gdac = tx.filter.H(fc/sim.fs);                    
       
%% Define how to calculate dc bias
% (i) full dc bias: sigtx*rclip
% (ii) reduced dc bias: sigtxdac*rclip
if sim.full_dc
    dc_bias = @(Pn) sim.rcliptx*sqrt(sum(2*Pn));
else
    dc_bias = @(Pn) sim.rcliptx*sqrt(sum(2*Pn.*abs(Gdac(nfc)).^2));
end

%% Power allocation
ofdm.power_allocation(tx, fiber, rx, sim);

% Signal std at transmitter and receiver
sigtx = sqrt(2*sum(ofdm.Pn));

% Check if power is at reasonable levels
assert(sigtx < 1e-3, 'Power is too high: %.2G (mW)\nCannot improve performance by increasing power.\n', 1e3*sqrt(2*sum(Pn)))

for navg = 1:sim.Navg

    xt = ofdm.generate_signal(tx, sim);

    %% Add dc bias (dc_bias determines whether it's added full dc bias or clipped dc bias)
    xt = xt + dc_bias(Pn);

    % Ensure that signal is non-negative
    xtc = xt;
    xtc(xt < 0) = 0;

    %% Apply frequency response of the laser
    [~, Et] = modulator(xtc, tx, sim);
           
    %% Fiber propagation
    Et = fiber.linear_propagation(Et);
    
    %% Direct detection
    Its = rx.R*abs(Et).^2;

    %% Add thermal noise
    wt = sqrt(rx.Sth*ofdm.fsct)*randn(1, sim.Mct*sim.Nsymb*(Npre_os + Nc));
    It = Its + wt;

    %% Add shot noise
    if sim.shot
        q = 1.60217657e-19;      % electron charge (C)
        Id = 0;                  % dark current

        % Instataneous received power considering only attenuation from the fiber   
        Sshot = 2*q*(rx.R*Ptz + Id);     % one-sided shot noise PSD

        % Frequency is divided by two because PSD is one-sided
        ws = sqrt(Sshot*ofdm.fsct/2).*randn(size(It));
        
        It = It + ws;
    end
    
    ber_temp = ofdm.detect(It);
    
    ber.count = [ber.count; ber_temp.count];
    ber.interval = [ber.interval; ber_temp.interval];
end

ber.count = mean(ber.count);
% ber.est = mean(ber.est);
ber.interval = mean(ber.interval, 2);



