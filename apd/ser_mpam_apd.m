clear, clc, close all

addpath f

M = [4 8 16];
colors = {'k', 'b', 'r'};

% M-PAM
mpam.Rb = 100e9;

% Simulation
sim.BERtarget = 1.8e-4;
sim.shot = 'on';
sim.rin = 'on';

simUniform = sim;
simUniform.levelSpacing = 'uniform';

simNonUniform = sim;
simNonUniform.levelSpacing = 'nonuniform';

% TX
PrecdBm = -30:0.1:0;
Prec = 1e-3*10.^(PrecdBm/10);
tx.RIN = -150;

% RX
% apd
apd.N0 = (30e-12)^2;
apd.R = 1;
apd.ka = 0.09;
apd.Gbw = 300e9;
apd.Fa = @(ka, G) ka*G + (1 - ka)*(2 - 1/G);

% pin
pin = apd;
pin.Gapd = 1;
pin.ka = 0;  

for k = 1:length(M)
    % M-PAM
    mpam.M = M(k);
    mpam.bw = mpam.Rb/log2(mpam.M);

    %% Uniform spacing
    % Gain constrained by the Bandwidth x Gain Product
    uapd_c = apd;
    uapd_c.Gapd = calcOptAPDGain(mpam, tx, uapd_c, simUniform, uapd_c.Gbw/mpam.bw);

    % Gain Unconstrained by the Bandwidth x Gain Product
    uapd_uc = apd;
    uapd_uc.Gapd = calcOptAPDGain(mpam, tx, uapd_uc, simUniform, 1000);

    % Calculate BER
    for kk = 1:length(Prec)
        tx.Prec = Prec(kk);
        % Constrained gain APD
        ber_u(k).apd_c(kk) = ber_mpam(mpam, tx, uapd_c, simUniform); 
        % Unconstrained gain APD
        ber_u(k).apd_uc(kk) = ber_mpam(mpam, tx, uapd_uc, simUniform); 
        % pin
        ber_u(k).pin(kk) = ber_mpam(mpam, tx, pin, simUniform); 
    end

    %% Non-uniform level spacing
    % Gain constrained by the Bandwidth x Gain Product
    napd_c = apd;
    napd_c.Gapd = calcOptAPDGain(mpam, tx, napd_c, simNonUniform, napd_c.Gbw/mpam.bw)

    [a, b] = calcOptLevelSpacing(mpam, tx, napd_c, simNonUniform);
    Popt_apd = mean(a);
    
    a/a(2);

    % Gain Unconstrained by the Bandwidth x Gain Product
    napd_uc = apd;
    napd_uc.Gapd = calcOptAPDGain(mpam, tx, napd_uc, simNonUniform, 1000)

    [a, b] = calcOptLevelSpacing(mpam, tx, napd_uc, simNonUniform);
    Popt_apd = mean(a);

    % pin
    [a, b] = calcOptLevelSpacing(mpam, tx, pin, simNonUniform);
    Popt_pin = mean(a);

    % Calculate BER
    for kk = 1:length(Prec)
        tx.Prec = Prec(kk);
        % Constrained gain APD
        ber_n(k).apd_c(kk) = ber_mpam(mpam, tx, uapd_c, simNonUniform); 
        % Unconstrained gain APD
        ber_n(k).apd_uc(kk) = ber_mpam(mpam, tx, uapd_uc, simNonUniform); 
        % pin
        ber_n(k).pin(kk) = ber_mpam(mpam, tx, pin, simNonUniform); 
    end


    %% Figures
    figure(1), hold on, box on, grid on
    plot(PrecdBm, log10(ber_u(k).apd_c),  '-', 'Color', colors{k}, 'LineWidth', 1.5)
    plot(PrecdBm, log10(ber_u(k).apd_uc),  ':', 'Color', colors{k}, 'LineWidth', 1.5)
    plot(PrecdBm, log10(ber_u(k).pin),  '--', 'Color', colors{k}, 'LineWidth', 1.5)
    xlabel('P_{rec} (dBm)', 'FontSize', 12)
    ylabel('log_{10}(BER)', 'FontSize', 12)
    legend('apd (maximum gain)', 'apd (optimal gain)', 'pin','Location', 'SouthWest')
    axis([-30 0 -7 0])

    figure(2), hold on, box on, grid on
    plot(PrecdBm, log10(ber_n(k).apd_c),  '-', 'Color', colors{k}, 'LineWidth', 1.5)
    plot(PrecdBm, log10(ber_n(k).apd_uc),  ':', 'Color', colors{k}, 'LineWidth', 1.5)
    plot(PrecdBm, log10(ber_n(k).pin),  '--', 'Color', colors{k}, 'LineWidth', 1.5)
    xlabel('P_{rec} (dBm)', 'FontSize', 12)
    ylabel('log_{10}(BER)', 'FontSize', 12)
    legend('apd (maximum gain)', 'apd (optimal gain)', 'pin','Location', 'SouthWest')
    axis([-30 0 -7 0])
end


