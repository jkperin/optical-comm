%% Simulation of analog coherent detection system
clear, clc

addpath ../f/
addpath ../analog/
addpath ../../f/

%% ======================== Simulation parameters =========================
sim.Nsymb = 2^14; % Number of symbols in montecarlo simulation
sim.Mct = 12;    % Oversampling ratio to simulate continuous time 
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 256; % number of symbols to be discarded from the begining and end of the sequence 
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.Rb = 2*112e9; % Bit rate
sim.Npol = 2;                                                              % number of polarizations
sim.pulse_shape = select_pulse_shape('rect', sim.Mct);                     % pulse shape
sim.ModFormat = QAM(4, sim.Rb/sim.Npol, sim.pulse_shape);                  % M-QAM modulation format
sim.Rs = sim.ModFormat.Rs;

sim.fs = sim.ModFormat.Rs*sim.Mct;
[f, t] = freq_time(sim.N, sim.fs);
dt = t(2) - t(1);

Analog.csi = 1/sqrt(2);                                                    % damping coefficient of second-order loop filter
Analog.Delay = 300e-12;                                                        % Additional loop delay in s (not including group delay from filters)

% Optimize EPLL parameters
Analog.wn = optimizePLL(Analog.csi, 1, Analog.Delay, 2*200e3, sim, false);
Analog.EPLL.nums = [2*Analog.csi*Analog.wn Analog.wn^2];
Analog.EPLL.dens = [1 0 0]; % descending powers of s
[Analog.EPLL.numz, Analog.EPLL.denz] = impinvar(Analog.EPLL.nums, Analog.EPLL.dens, sim.fs);
fprintf('Loop filter fn: %.3f GHz\n', Analog.wn/(2*pi*1e9));

% Loop filter
LoopFilter = ClassFilter(Analog.EPLL.numz, Analog.EPLL.denz, sim.fs);

Tstep = sim.N/4; % instant where step in frequency occurs
Fstep = (0.1:0.1:1)*1e9;
tol = 1e-3;
for n = 1:length(Fstep)
    fstep = Fstep(n);
    freqOffset = [zeros(1, Tstep-1) fstep*ones(1, sim.N-Tstep+1)];                                                   % Frequency shift with respect to transmitter laser in Hz
    phiS = 2*pi*freqOffset.*t;
    y = zeros(size(phiS));
    phiVCO = zeros(size(phiS));
    LoopFilter = LoopFilter.reset();
    for k = 2:length(t)
%        y(k) = -(sin(phiS(k) - phiVCO(k-1))); % Costas
       y(k) = -sign(sin(phiS(k) - phiVCO(k-1))); % Logic

       phiVCO(k) = LoopFilter.filter(y(k));
    end
    
    % Wrapped difference between signal and VCO phase
    phaseError = phiS-phiVCO;
    phaseErrorf = filtfilt(ones(1, 101)/101, 1, phaseError);
    dphi = asin(sin(phaseError));

    % Measure time until locking again
    adphi = abs(diff(phaseErrorf));
    S = cumsum(adphi);
    S = S/S(end);
    plock = find(S >= 1-tol, 1, 'first');
    tlock(n) = t(plock) - t(Tstep);
    
    % Count cycle slips
    [peaks, idx] = findpeaks(dphi);
    Ncycle_slips(n) = sum(peaks >= 0.9*pi/2);
    idx = idx(peaks >= 0.9*pi/2);
    
    fVCO = diff(phiVCO)/(2*pi*dt*1e6);
    figure(1), clf
    subplot(211), box on, hold on
    plot(t*1e9, freqOffset/1e6, 'k')
    plot(t(1:end-1)*1e9, fVCO)
    plot(t(idx)*1e9, fVCO(idx), 'xk')
    a = axis;
    plot(1e9*t(plock)*[1 1],  a(3:4), ':k')
    xlabel('Time (ns)')
    ylabel('Frequency (MHz)')
    axis tight
    subplot(212), hold on, box on
    plot(t*1e9, dphi)
    plot(t(idx)*1e9, dphi(idx), 'xk')
    a = axis;
    plot(1e9*t(plock)*[1 1],  a(3:4), ':k')
    xlabel('Time (ns)')
    ylabel('Wrapped phase error (rad)')
    drawnow
end

Colors = {[51, 105, 232]/255; [228,26,28]/255; [0,153,37]/255; [255,127,0]/255};

if length(Fstep) > 1
    figure(2), box on, hold on
%     child = allchild(gca);
    plot(Fstep/1e6, tlock*1e9, '--s', 'Color', Colors{3})
    ylabel('Recovery time (ns)', 'FontSize', 12)
    xlabel('Frequency step (MHz)', 'FontSize', 12)
    set(gca, 'FontSize', 12)

    figure(3), box on, hold on
    plot(Fstep/1e6, Ncycle_slips, '--s', 'Color', Colors{3})
    xlabel('Frequency step (MHz)')
    ylabel('Number of cycle slips')
end


    