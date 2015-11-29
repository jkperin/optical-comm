%% Calculate noise distribution using doubly-stochastic random process, gaussian approximation, and heuristic pdf
clear, clc, close all

format compact
addpath f

%% Simulation
Nsamples = 2^12;

%% M-PAM
mpam.M = 4;
mpam.Rb = 100e9;
mpam.Rs = mpam.Rb/log2(mpam.M);

%% APD
% class apd(GainBW, ka, Gain, noise_stats)
apd = apd(10*log10(300e9/mpam.Rs), 0.09, 300e9);

PlevelsdBm = -30;
Plevels = 1e-3*10.^(PlevelsdBm/10);

pdfs = apd.levels_pdf(Plevels, mpam.Rs);

for k = 1:length(Plevels)
   
    Pin = Plevels(k)*ones(Nsamples, 1);

    % Get heuristic pdfs
    out_gauss = apd.detect(Pin, mpam.Rs, 'gaussian');
    out_doubl = apd.detect(Pin, mpam.Rs, 'gaussian');
    
    out_doubl = out_doubl + sqrt(mpam.Rs/2*(20e-12)^2)*randn(size(out_doubl));

    [p_hist, Ihist] = hist(out_doubl, 50);
    p_hist = p_hist/trapz(Ihist, p_hist);
    
    [px, x] = apd.output_pdf_saddlepoint(Plevels(k), mpam.Rs, (20e-12)^2);
    for kk = 1:length(x)
        [cdf_left(kk), shat_left(kk)] = apd.tail_saddlepoint_approx(x(kk), Plevels(k), mpam.Rs, (20e-12)^2, 'left');
        [cdf_right(kk), shat_right(kk)] = apd.tail_saddlepoint_approx(x(kk), Plevels(k), mpam.Rs, (20e-12)^2, 'right');
    end
    
    figure(1), hold on
    plot(1e3*pdfs(k).I, pdfs(k).p)
    plot(1e3*pdfs(k).I, pdfs(k).p_gauss)
    bar(1e3*Ihist, p_hist)
    plot(x*1e3, px, 'k')
    legend('doubly-stochastic', 'gaussian', 'histogram')
    xlabel('Current (mA)', 'FontSize', 12)
    ylabel('pdf', 'FontSize', 12)
    axis([0 1e3*pdfs(k).I(end) 0 5e5])

    figure, hold on
    plot(x, cumtrapz(x, px))
    plot(x, cdf_left)
    plot(x, cdf_right)
    axis([x(1) x(end) 0 1])
    legend('true', 'left', 'right')
    
    figure, hold on
    plot(x, shat_left)
    plot(x, shat_right)

%     figure(2), hold on
%     plot(1e3*pdfs(k).I, log10(pdfs(k).p), 'k')
%     plot(1e3*pdfs(k).I, log10(pdfs(k).p_gauss), '--r')
%     legend('doubly-stochastic', 'gaussian', 'Location', 'SouthEast')
%     xlabel('Current (mA)', 'FontSize', 12)
%     ylabel('log(pdf)', 'FontSize', 12)
%     axis([0 0.08 -10 8])
end


