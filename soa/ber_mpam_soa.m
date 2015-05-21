clear, clc, close all

addpath f

M = [4 8 16];
colors = {'k', 'b', 'r'};

% Constants
h = 6.62606957e-34;
q = 1.60217657e-19;
c = 299792458;

% M-PAM
mpam.Rb = 100e9;

% Simulation
sim.BERtarget = 1.8e-4; % Target BER and SER
sim.Nsymb = 2^14;
sim.Mct = 1;    % Oversampling ratio for continuous-time simulation
sim.soa = true;
sim.shot = 'on';
sim.rin = 'on';
sim.levelSpacing = 'nonuniform';

% simUniform = sim;
% simUniform.levelSpacing = 'uniform';
% 
% simNonUniform = sim;
% simNonUniform.levelSpacing = 'nonuniform';

% TX
PrecdBm = -30:0;
Prec = 1e-3*10.^(PrecdBm/10);
tx.RIN = -150;
tx.lamb = 1310e-9;

% RX
% pin
rx.N0 = (30e-12)^2;
rx.R = 1;
rx.ka = 0;
rx.Gapd = 1;
rx.Fa = @(ka, G) ka*G + (1 - ka)*(2 - 1/G);

% APD
apd.N0 = (30e-12)^2;
apd.R = 1;
apd.ka = 0.09;
apd.Gbw = 300e9;
apd.Fa = @(ka, G) ka*G + (1 - ka)*(2 - 1/G);

% SOA
soa.Fn = 9; % dB
soa.Seq = @(G) (G - 1)*10^(soa.Fn/10)/2*(h*c/tx.lamb); % one-sided PSD
soa.G = 10^(5/10);


for k = 1:length(M)
    % M-PAM
    mpam.M = M(k);
    mpam.bw = mpam.Rb/log2(mpam.M);
    mpam.Rs = mpam.Rb/log2(mpam.M);

   % Generate unipolar PAM signal
    if strcmp(sim.levelSpacing, 'uniform')
        mpam.a = (0:2:2*(mpam.M-1)).';
        mpam.b = (1:2:(2*(mpam.M-1)-1)).';
    elseif strcmp(sim.levelSpacing, 'nonuniform')
        % Optimal level spacing for APD is calculated in ber_mpam.m
        [mpam.a, mpam.b] = calcOptLevelSpacingAmpl(mpam, tx, soa, rx, sim);   
    else
        error('Invalide Option!')
    end 

    % APD
    apd.Gapd = calcOptAPDGain(mpam, tx, apd, sim, 100)
%     apd.Gapd = calcOptAPDGain(mpam, tx, apd, sim, apd.Gbw/mpam.Rs);
    
    % SOA
    soa.G = calcOptSOAGain(mpam, tx, soa, rx, sim, 100)

    % Calculate BER
    for kk = 1:length(Prec)
        tx.Prec = Prec(kk);

        mpam.b = mpam.b*tx.Prec/mean(mpam.a);
        mpam.a = mpam.a*tx.Prec/mean(mpam.a);

        % theoretical BER
        berapd(kk) = ber_mpam(mpam, tx, apd, sim);

        bersoa(kk) =  mpam_ber_soa(mpam, soa, rx);
        
        % Optimal level spacing for APD is calculated in this function
        berpin(kk) = ber_mpam(mpam, tx, rx, sim);
    end

    %% Figures
    figure(1), hold on, box on, grid on
    plot(PrecdBm, log10(berapd), '-', 'Color', colors{k}, 'LineWidth', 1.5)
    plot(PrecdBm, log10(bersoa), ':', 'Color', colors{k}, 'LineWidth', 1.5)
    plot(PrecdBm, log10(berpin), '--', 'Color', colors{k}, 'LineWidth', 1.5)
    xlabel('P_{rec} (dBm)', 'FontSize', 12)
    ylabel('log_{10}(BER)', 'FontSize', 12)
    legend('APD', 'SOA', 'PIN', 'Location', 'SouthWest')
    axis([-30 0 -7 0])

%     figure(2), hold on, box on, grid on
%     plot(PrecdBm, log10(ber_n(k).apd_c),  '-', 'Color', colors{k}, 'LineWidth', 1.5)
%     plot(PrecdBm, log10(ber_n(k).apd_uc),  ':', 'Color', colors{k}, 'LineWidth', 1.5)
%     plot(PrecdBm, log10(ber_n(k).pin),  '--', 'Color', colors{k}, 'LineWidth', 1.5)
%     xlabel('P_{rec} (dBm)', 'FontSize', 12)
%     ylabel('log_{10}(BER)', 'FontSize', 12)
%     legend('apd (maximum gain)', 'apd (optimal gain)', 'pin','Location', 'SouthWest')
%     axis([-30 0 -7 0])
end


