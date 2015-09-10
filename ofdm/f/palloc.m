function [Pn, CS] = palloc(ofdm, tx, fiber, rx, sim)
Nit_piterative = 50;            % Number of iterations used to estimate required power at subcarriers to achieve target BER for the case of quantization.

Nc = ofdm.Nc; 
Nu = ofdm.Nu;
fc = ofdm.fc;
CS = ofdm.CS*ones(1,Nu/2);

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

if isfield(sim, 'full_dc') && sim.full_dc
    dc_bias = @(Pn) tx.rclip*sqrt(sum(2*Pn));
else
    dc_bias = @(Pn) tx.rclip*sqrt(sum(2*Pn.*abs(Gdac).^2));
end

% Iterate to get correct power. For kk = 1, it assumes that the
% power is zero i.e., there is no noise due to quantization
varQ1th = 0;
varQ2th = 0;
varshot = 0;
varrin = 0;
kk = 1;
Pchange = 1;
Pn1 = 0;
% Iterante until power changes is less than .1% or when maximum
% number of iterations is reached.
while Pchange > 1e-3 && kk < Nit_piterative
    % Gain-to-noise ratio
    GNR = Nc*abs(Gch).^2./(ofdm.fs*rx.Sth*abs(Gadc).^2 + varQ1th*abs(Gch).^2 + varQ2th + varshot + varrin); 

    % Run Levin-Campello algorithm to determine optimal power allocation
    [~, CS, Pn] = Levin_Campello_MA(ofdm.B, 1, sim.BERtarget, GNR); 

    % Signal std at transmitter and receiver
    sigtx = sqrt(2*sum(Pn));
    sigrx = sqrt(2*sum(Pn.*abs(Gch).^2));
    Ptx_est = dc_bias(Pn);

    %% Add additional dc-bias due to finite exctinction ratio
    rex = 10^(tx.rexdB/10);

    Pmin_rex = 2*tx.rclip*sigtx/(rex - 1); %%additional dc bias

    Ptx_est = Ptx_est + Pmin_rex; 
    
    %% Quantization 
    if isfield(sim, 'quantiz') && sim.quantiz
        % Quantization noise variance
        delta1th = 2*tx.rclip*sigtx/(2^sim.ENOB-1);
        delta2th = 2*rx.rclip*sigrx/(2^sim.ENOB-1); 
        varQ1th = (1 - 2*qfunc(tx.rclip))*delta1th^2/12; % quantization at the transmitter
        varQ2th = (1 - 2*qfunc(rx.rclip))*delta2th^2/12; % quantization at the receiver
    end

    %% Shot noise
    if isfield(sim, 'shot') && sim.shot
        q = 1.60217657e-19;      % electron charge (C)
        Id = 0;                  % dark current

        Prx_est = Ptx_est*fiber.link_attenuation(tx.lamb);

        % Shot noise psd (one-sided)
        Sshot = 2*q*(rx.R*Prx_est + Id);

        % Shot noise variance. ofdm.fs/2 because Sshot is one-sided
        varshot = Sshot*ofdm.fs/2*abs(Gadc).^2;    
    end

    %% RIN
    if isfield(sim, 'RIN') && sim.RIN         
        % Calculate double-sided RIN PSD, which depends on the average power
        Srin = 10^(tx.RIN/10)*Ptx_est.^2;

        % RIN variance. sim.fs because Srin is double-sided.
        varrin = Srin*ofdm.fs.*abs(Hfiber.*rx.R.*Gadc).^2;
    end  

    Pchange = sum(abs(Pn - Pn1)./Pn);

    Pn1 = Pn;

    kk = kk + 1;
end

% Auxiliary variables. They're used to check approximations and validate
% the code
if isfield(sim, 'verbose') && sim.verbose && isfield(sim, 'quantiz') && sim.quantiz
    aux.delta1th = delta1th;
    aux.delta2th = delta2th;
    aux.varQ1th = varQ1th;
    aux.varQ2th = varQ2th;
    
    aux
end