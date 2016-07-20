function [ysamp, idx] = analog_time_recovery(yct, TimeRec, sim, verbose)
%% Time recovery for analog-based receiver
% Inputs:
% - yc: input signal in "continuous-time". If yct is complex, ony the I 
% component will be used for estimating the clock.
% - TimeRec: struct containing time recovery method properties
%     - type: eithehr 'squarer' or 'Mueller-Muller'
%     - bpf (only if type = 'squarer'): band pass filter struct 
%     - csi (type = 'Mueller-Muller'): damping of loop filter
%     - wn (type = 'Mueller-Muller'): relaxation frequency of loop filter
%     - CT2DT (type = 'Mueller-Muller', optional, default = 'bilinear')
%     method of converting continuous-time loop filter into discrete time
% - sim: struct containing simulation parameters
% Outputs:
% - ysamp: sampled signal
% - idx: indices such that ysamp = yct(idx)

switch lower(TimeRec.type)
    case 'spectral-line-bpf'
        %% Spectral line syncronizer based on squaring and bandpass filtering
        % Differentiator -> nonlinearity (squaring) -> band pass filter
        % centered at the symbol rate
        % Note: BPF could be replaced by a PLL, but this is not considered
        % here
        
        yct = resample(yct, TimeRec.Mct, sim.Mct);
        fs = sim.fs*TimeRec.Mct/sim.Mct;
        [f, t] = freq_time(length(yct), fs);
        Hbpf = ifftshift(TimeRec.bpf.H(f/fs));
        
        % Nonlinearity + differentiator
        y = real(yct);
%         ynl = diff([y 0]).^2;
        ynl = y.^2;
        
        % BPF
        clk = real(ifft(fft(ynl).*Hbpf));
        
        % Limiting amps
        clk(clk > 0) = 1;
        clk(clk < 0) = 0;
        
        % Rising-edge detection
        samp = diff([clk 0]);
        idx = find(samp > 0);
        
        % Sampling
        tsamp = t(idx);
        ysamp = trim_pad(yct(idx));        
        
        % Plots
        if exist('verbose', 'var') && verbose
            figure(405), clf, hold on
            plot(1e12*(sim.t(1:sim.Mct:min(length(sim.t), length(tsamp)*sim.Mct)) - tsamp))
            xlabel('Symbol')
            ylabel('Sampling time error (ps)')
            title('Timing recovery')
            scatterplot(ysamp)
        end
 
    case 'spectral-line-pll'
        ysamp = yct(1:sim.Mct:end); 
    case 'none'
        ysamp = yct(1:sim.Mct:end);     
    otherwise
        error('analog_time_recovery: invalide time recovery type')
end

function y = trim_pad(x)
%% Ensures that vector x has sim.Nsymb points by trimming or padding zeros
    N = sim.Nsymb;
    if length(x) > N
        y = x(1:N);
    elseif length(x) < N
        y = [x zeros(1, N-length(x))];
    else
        y = x;
    end
end
end
