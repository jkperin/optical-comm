%% Maximum delay for a given receiver sensitivity penalty as a function of total linewidth
% Loop filter bandwidth is optimized for each linewidth and delay
clear, clc

addpath ../../f/
addpath ../f/

Qpsk = QAM(4, 2*112e9); % Modulation format
BERtarget = 1.8e-4; % target BER
sim.BERtarget = BERtarget;
sim.ModFormat = Qpsk;
SNRdBref = SNRreq(BERtarget, Qpsk.M, 'QAM'); % Required SNR to achieve target BER
R = 1; % photodiodes responsivity
q = 1.60217662e-19; % electron charge
csi = sqrt(2)/2; % damping factor or loop filter
wnOffFactor = 1; % factor by which optimal loop filter bandwidth is multiplied

totalLinewidth = 2*1e3*(100:25:1e3);

% Power penalty assuming ideal Prx as reference
Ppen = 0.5;

for Ncpr = 1:2 % number of polarizations used in CPR 
    for k = 1:length(totalLinewidth)
        [maxDelay(Ncpr, k), wnOpt(Ncpr, k)] = max_loop_delay(SNRdBref + Ppen, csi, totalLinewidth(k), Ncpr, sim, wnOffFactor);
    end
end

%
Line = {'--', '-'};
figure(3), hold on, box on
plot(totalLinewidth/1e3, 1e12*maxDelay(1, :), Line{1}, 'LineWidth', 2, 'DisplayName', ['1 pol, \omega_n = ' num2str(wnOffFactor) '\omega_{n, opt}'])
plot(totalLinewidth/1e3, 1e12*maxDelay(2, :), Line{2}, 'LineWidth', 2, 'DisplayName', ['2 pols, \omega_n = ' num2str(wnOffFactor) '\omega_{n, opt}'])
xlabel('Total linewidth (kHz)', 'FontSize', 12)
ylabel('Maximum delay for 0.5-dB power penalty (ps)', 'FontSize', 12)
set(gca, 'FontSize', 12)
legend('-DynamicLegend')
axis([totalLinewidth([1 end])/1e3 0 2000])
grid on
% m2tikz = matlab2tikz(gca);
% m2tikz.write('max_delay_vs_linewidth.tikz');

figure(4), clf, hold on, box on
plot(totalLinewidth/1e3, wnOpt(1, :)/(2*pi*1e6), Line{1}, 'LineWidth', 2)
plot(totalLinewidth/1e3, wnOpt(2, :)/(2*pi*1e6), Line{2}, 'LineWidth', 2)
xlabel('Total linewidth (kHz)', 'FontSize', 12)
ylabel('Optimal loop filter bandwidth (MHz)', 'FontSize', 12)
set(gca, 'FontSize', 12)
axis([totalLinewidth([1 end])/1e3 0 200])
legend('1 pol', '2 pols')
grid on
% m2tikz2 = matlab2tikz(gca);
% m2tikz2.write('omegan_opt_vs_linewidth.tikz');

