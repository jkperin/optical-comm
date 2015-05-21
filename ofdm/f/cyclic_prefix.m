%% Calculate the cyclic prefix length after oversampling 
% The channel is assumed to be a 2nd-order filter with damping ratio = 1.

% Algorithm:
% 1. Assumes a cyclic prefix length k
% 2. From k, it calculates the new sampling rate fs
% 3. Calculate the number of samples at both sides (Nneg and Npos) that 
% contains the desired fraction of energy
% 4. From that calculates Npre_os = Npos + Neg; (zero is not included)
% 5. If Npre_os == k end simulation, otherwise increment k and repeat

function [Npre_os, Nneg, Npos] = cyclic_prefix(ofdm, tx, rx, sim)

frac_incl = 0.9999;    % fraction of energy to be included within CP length
max_it = 100;          % maximum number of iterations
Mct = sim.Mct;
Ntot = 1024;

switch ofdm.ofdm
    case 'aco_ofdm'
        Nd = 2*ofdm.Nu;  
    case 'dc_ofdm'
        Nd = ofdm.Nu;
    otherwise
        error('Invalid option!')
end

k = 0;
Npre_os = Inf;
while Npre_os > k && k < max_it
    fs = ofdm.Rs*(ofdm.Nc + k)/Nd;
    fsct = Mct*fs;
    dt = 1/fsct;
    
    % Group delay of modulator in samples
    hl_delay = tx.hl_delay*fsct;
    
    % Channel impulse response
    tct = (0:Ntot-1)*dt;
    hct = tx.hl(tct);
    gdac = tx.gdac;
    gadc = rx.gadc;
    
    % Total impulse response in continuous time
    pct = conv(conv(gdac, hct, 'full'), gadc, 'full');
    pct = pct/max(pct);

    % Remove group delay due to filters so that impulse response is
    % centered at zero.
    tct = (0:length(pct)-1)*dt;
    tct = tct - ceil(tx.gdac_delay + rx.gadc_delay + hl_delay)*dt;
    
    n0 = ceil(tx.gdac_delay + hl_delay + rx.gadc_delay) + 1; % new zero after remove group delay
    
    % Sampling at the chip rate
    t = tct([fliplr(n0:-Mct:1), n0+Mct:Mct:end]);
    p = pct([fliplr(n0:-Mct:1), n0+Mct:Mct:end]);
       
    % CP based on energy
    en_frac = cumsum(p.^2)/sum(p.^2);
    Nneg = sum(en_frac >= (1 - frac_incl)/2 & t < 0);
    Npos = sum(en_frac <= (1 + frac_incl)/2 & t > 0);
    
    % Number of samples necessary to attain desired fraction of energy
    Npre_os = Nneg + Npos;

    k = k + 1;
end


%% plot and verify if cyclic prefix contains desired portion of energy
% figure
% plot(tct, pct)
% hold on
% stem(t, p, 'fill')
% plot(Mct*Npos*[1 1]/fsct, [-1 1], 'k')
% plot(-Mct*Nneg*[1 1]/fsct, [-1 1], 'k')
% axis([-Mct*Npre_os/fsct Mct*Npre_os/fsct -0.5 1])
% title(sprintf('Cutoff Frequency: %d', tx.fnl/1e9))
% 

assert(k ~= max_it, 'CP calculation did not converge');