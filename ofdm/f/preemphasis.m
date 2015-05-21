function [Pn, CS, aux] = preemphasis(ofdm, tx, fiber, rx, Gch, sim)
Nit_piterative = 50;            % Number of iterations used to estimate required power at subcarriers to achieve target BER for the case of quantization.

Nc = ofdm.Nc; 
Nu = ofdm.Nu;

CS = ofdm.CS*ones(1,Nu/2);

K = 1 - 2*qfunc(sim.rcliptx);                % Amplitude attenuation due to clipping = (1-Q(r))

% SNR to achieve target BER sim.Pb
snrdB = fzero(@(x) berqam(ofdm.CS, x) - sim.Pb, 20);
snrl = 10^(snrdB/10);

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
    Pnrx = snrl*(ofdm.fs*rx.Sth*abs(rx.Gadc).^2 + varQ1th*abs(Gch).^2 + varQ2th + varshot + varrin)/Nc;                                               

    Pn = Pnrx./(K^2*abs(Gch).^2);

    % Signal std at transmitter and receiver
    sigtx = sqrt(2*sum(Pn));
    sigrx = sqrt(2*sum(Pn.*K^2.*abs(Gch).^2));
    Ptx_est = tx.kappa*ofdm.dc_bias(Pn);
    
    %% Add additional dc-bias due to finite exctinction ratio
    if isfield(tx, 'rexdB') % check if exctinction ratio (rex) was defined
        rex = 10^(tx.rexdB/10);
        
        Pmin_rex = 2*sim.rcliptx*sigtx/(rex - 1); %% additional dc bias
        
        Ptx_est = Ptx_est + Pmin_rex; 
                                                
    end
    

    if sim.quantiz
        % Quantization noise variance
        delta1th = 2*sim.rcliptx*sigtx/(2^sim.ENOB-1);
        delta2th = 2*sim.rcliprx*sigrx/(2^sim.ENOB-1); 
        varQ1th = (1 - 2*qfunc(sim.rcliptx))*delta1th^2/12; % quantization at the transmitter
        varQ2th = (1 - 2*qfunc(sim.rcliprx))*delta2th^2/12; % quantization at the receiver
    end

    if sim.shot
        q = 1.60217657e-19;      % electron charge (C)
        Id = 0;                  % dark current

        Prx_est = Ptx_est/10^(fiber.att*fiber.L/1e4);

        % Shot noise psd (one-sided)
        Sshot = 2*q*(rx.R*Prx_est + Id);

        % Shot noise variance. ofdm.fs/2 because Sshot is one-sided
        varshot = Sshot*ofdm.fs/2*abs(rx.Gadc).^2;    
    end

    if sim.RIN                
        % Calculate double-sided RIN PSD, which depends on the average power
        Srin = 10^(tx.RIN/10)*Ptx_est.^2;

        % RIN variance. ofdm.fs because Srin is double-sided.
        varrin = Srin*ofdm.fs.*abs(10^(-fiber.att*fiber.L/1e4)*fiber.Hcd.*rx.R.*rx.Gadc).^2;
    end                

    Pchange = sum(abs(Pn - Pn1)./Pn);

    Pn1 = Pn;

    kk = kk + 1;
end  

% Auxiliary variables. They're used to check approximations and validate
% the code
aux.delta1th = delta1th;
aux.delta2th = delta2th;
aux.varQ1th = varQ1th;
aux.varQ2th = varQ2th;