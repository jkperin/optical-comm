function [SignalOut, ASEf, GaindB, GFFdB] = Amp_GFF_Loss(E, Pump, Signal, ASEf, lossdB, OptPowerProfiledBm, GFFdB)
%% Propagete Signal over amplifier, GFF, and some loss 
% Inputs:
%   - E: instance of class EDF
%   - Pump: instance of class Channels corresponding to the Pump
%   - Signal: instance of class Channels corresponding to the input singal
%   to the amplifier
%   - ASEf: instance of class Channels corresponding to the forward ASE
%   - lossdB: loss profile in dB. Loss could be due to a fiber span, for
%   instance.
%   - OptPowerProfiledBm: optimized power profile
%   - GFFdB (optional, default = optimal): gain profile of gain flatening
%   filter. If not provided, GFFdB is the optimal GFF.
% Outputs:
%   - SignalOut: instance of class Channels corresponding to the signals at
%   the output of fiber, i.e., after amp, GFF, and fiber.
%   - ASEf: instance of class Channels corresponding to the foward ASE at
%   the output of the fiber.
%   - GaindB: amplifier gain profile in dB. This is computed using the
%   two-level model
%   - GFFdB: gain profile of gain flatening filter.

    % Constants
    OffPower = 1e-12; % power set for OFF channels

    SignalOut = Signal;
    
    %% Amplifier
    ASEb = Channels(Signal.wavelength, 0, 'backward'); % set backward ASE to zero at every iteration

    [GaindB, ~, Psignal_out, Pase, sol] = E.propagate(Pump, Signal, ASEf, ASEb, ASEf.df, 'two-level', 50, false);

    if any(sol.y(:) < 0)
        warning('EDF/propagate: solution contains negative power')
    end

    PoutdBm = Watt2dBm(Psignal_out);
    
    %% GFF and span attenuation
    if not(exist('GFF', 'var')) % use optimal GFF filter
        % GFF gain profile is negative dB, whereas spanAttdB is positive dB
        GFFdB = min(OptPowerProfiledBm - PoutdBm + lossdB, 0);
    end
    
    SignalOut.PdBm = PoutdBm + GFFdB - lossdB;
    SignalOut.P(SignalOut.P < OffPower) = OffPower; % Set power level to minimum OffPower to ensure numerical stability 
    
    ASEf.PdBm = Watt2dBm(Pase) + GFFdB - lossdB;

end