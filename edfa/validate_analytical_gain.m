clear, clc

% Pump parameters
Pump.P = 5e-3;
Pump.wavelength = 1480e-9;
Pump.alpha = 10^(1.6/10);
Pump.gain_coeff = 10^(0.5/10);
Pump.direction = 1;

% Signal parameters
Signal.P = [100e-9 100e-9];
Signal.wavelength = [1545e-9 1550e-9];
Signal.alpha = [10^(2.6/10) 10^(2.6/10)];
Signal.gain_coeff = [10^(3.6/10) 10^(3.6/10)];
Signal.direction = [1 1];

Amp = EDFA(Pump, Signal, 20);

Ppump = linspace(1, 10, 20)*1e-3;
Gain = zeros(length(Ppump), length(Signal.wavelength));
for k = 1:length(Ppump) 
    Amp.Pump.P = Ppump(k);
    Gain(k, :) = Amp.analytical_gain();
end

figure(1), box on, hold on
plot(Ppump*1e3, 10*log10(Gain))
xlabel('Pump power (mW)', 'FontSize', 12)
ylabel('Gain (dB)', 'FontSize', 12)
 axis([0 10 -10 40])
 
 Amp.optimal_length()
 
[ASE, sol] = Amp.amplified_power();

sol

figure, hold on, box on
l = ASE.wavelength(1:length(ASE.wavelength)/2);
plot(l*1e9, ASE.PoutdBm_b(:, end))
plot(l*1e9, ASE.PoutdBm_f(:, end))