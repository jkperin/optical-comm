function Y = dual_pol_stokes_receiver(Erec, Rx, sim)
%% Dual polarization coherent receiver: includes 90deg hybrid, 4 balanced photodiodes, transimpedance amplifier adding thermal noise, and coarse automatic gain control
% Inputs:
% - Erec: electric field of received signal [2 x N]
% - ELO: electric field of local oscillator [2 x N]
% - Rx: receiver parameters
% - sim: simulation parameters
% Note: if Erec and ELO are for [2 x 1] i.e., OPLL simulation, then this
% function and subfunctions cannot performed filtering

% Power splitters before 90-degree hybrid
[Es1, EX1] = powerSplitter(Erec(1, :), 0, 0.5, 0);
[Es2, EY1] = powerSplitter(Erec(2, :), 0, 0.5, 0);

% 90deg Hybrid
% Note: fields are 2 x N matrices, where each row is one polarization
% Power splitter of signal and LO for pol 1
[Es1i, Es1q] = powerSplitter(Es1, 0, 0.5, 0);
% Power splitter of signal and LO for pol 2
[Es2i, Es2q] = powerSplitter(Es2, 0, 0.5, 0);

% Note: last indices 1 and 2 refer to output ports rather than polarizations
[E1i1, E1i2] = powerSplitter(Es1i, Es2i*exp(-1j*pi/2), 0.5, 0);
[E1q1, E1q2] = powerSplitter(Es1q, Es2q, 0.5, 0);

% perform photodetection
Y1 = Rx.PD.detect(EX1, sim.fs, 'gaussian'); % |Erec(1, :)|^2

Y2 = Rx.PD.detect(EY1, sim.fs, 'gaussian'); % |Erec(2, :)|^2

Y3 = Rx.PD.detect(E1i1, sim.fs, 'gaussian')...
    - Rx.PD.detect(E1i2, sim.fs, 'gaussian'); % Re{Erec(1, :).*conj(Erec(2, :))}

Y4 = Rx.PD.detect(E1q1, sim.fs, 'gaussian')...
    - Rx.PD.detect(E1q2, sim.fs, 'gaussian'); % Im{Erec(1, :).*conj(Erec(2, :))}

% Add thermal noise
Y1 = Y1 + sqrt(Rx.N0*sim.fs/2)*randn(size(Y1));
Y2 = Y2 + sqrt(Rx.N0*sim.fs/2)*randn(size(Y2));
Y3 = Y3 + sqrt(Rx.N0*sim.fs/2)*randn(size(Y3));
Y4 = Y4 + sqrt(Rx.N0*sim.fs/2)*randn(size(Y4));

% Form output
Y = [Y1; Y2; Y3; Y4];