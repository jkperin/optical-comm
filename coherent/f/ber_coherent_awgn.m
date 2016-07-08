function [ber_ideal, ber_noise_enhancement, SNRdBtheory] = ber_coherent_awgn(Tx, Fiber, Rx, sim)
%% Estimate ideal BER of a coherent detection system
% Inputs:
% - Tx: transmitter parameters {filt: transmiter filter, Mod: modulator
% parameters}
% - Fiber: fiber parameters
% - Rx: receiver parameters
% - sim: simulation parameters {M: modulation order, ModFormat: modulation
% format, Rs: symbol rate, Mct: oversampling ratio to emulate continuous
% time, f: frequency vector}
% Outputs:
% - ber_ideal: BER of a system whose receiver filter is matched to the
% transmitter pulse shape. This doesn't take into account any frequency
% response of the channel and noise enhancement penalty
% - ber_noise_enhancement: BER of a system whose receiver filter is matched 
% to the receiver pulse shape and there is symbol rate liner equalization. 

% BER in AWGN channel
berAWGN = @(SNRdB) berawgn(SNRdB - 10*log10(2) - 10*log10(log2(sim.M)), lower(sim.ModFormat), sim.M);
% Note: - 10*log10(2) is due to the fact SNR = 2Es/N0

% Power levels
Prx = dBm2Watt(Tx.PlaunchdBm)/Fiber.link_attenuation(Tx.Laser.lambda); % received power
Plo = dBm2Watt(Rx.LO.PdBm);                                            % local oscillator power
Ppd = abs(sqrt(Plo/(4*sim.Npol)) + sqrt(Prx/(4*sim.Npol))).^2;         % incident power in each photodiode
Psig = Plo*Prx/(2*sim.Npol*sim.Npol);                          % Signal power per real dimension

% Estimate SNR of an ideal symbol-rate linear equalizer   
noiseBW = sim.Rs/2;                                                    % noise bandwidth of matched filter. This doesn't include noise enhancement penalty
varShot = 2*Rx.PD.varShot(Ppd, noiseBW);                               % Shot noise variance per real dimension
varThermal = Rx.N0*noiseBW;                                            % Thermal noise variance per real dimension
SNRdBtheory = 10*log10(2*Psig./(varShot + varThermal));
ber_ideal = berAWGN(SNRdBtheory);

% Design linear equalizer to estimate noise enhancement penalty
% Calculate noise bandwidth of an ideal receiver consisting of a matched 
% filter and symbol-rate linear equalization. This is only used to
% estiamate theoretical BER
eq.type = 'fixed td-sr-le';
eq.Ntaps = 31;
pulse_shape = select_pulse_shape('rect', sim.Mct);
mpam = PAM(sqrt(sim.M), sim.Rs, 'equally-spaced', pulse_shape);
Hch = (Tx.filt.H(sim.f/sim.fs).*Tx.Mod.Hel.*Fiber.Hdisp(sim.f, Tx.Laser.lambda)).';
sim.f = sim.f.';
[~, eq] = equalize(eq, [], Hch, mpam, sim); % this generates zero-forcing linear equalizer, which should be a good
% approximation of MMSE linear equalizer.
noiseBW_noiseEnhance = trapz(sim.f, abs(eq.Hrx.*eq.Hff(sim.f/sim.Rs)).^2)/2;

% Estimate SNR including noise enhacement penalty of an ideal symbol-rate linear equalizer   
varShot = 2*Rx.PD.varShot(Ppd, noiseBW_noiseEnhance);                               % Shot noise variance per real dimension
varThermal = Rx.N0*noiseBW_noiseEnhance;                                            % Thermal noise variance per real dimension
SNRdB_noise_enhancement = 10*log10(2*Psig./(varShot + varThermal));
ber_noise_enhancement = berAWGN(SNRdB_noise_enhancement);