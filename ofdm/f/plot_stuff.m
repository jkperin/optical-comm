%% Plot more important things 
% this script is called by the functions dc_ofdm and aco_ofdm to plot more
% important things

% Anonymous functions
AmpNorm = @(x) x/max(abs(x));  % Normalize amplitude (used for plotting)


%% Impulse and frequency response
dt = 1/ofdm.fsct;
tct = 0:dt:4*Npre_os/ofdm.fs-dt;
hct = tx.hl(tct);

% Total impulse response in continuous time
pct = conv(conv(tx.gdac, hct, 'full'), rx.gadc, 'full');
pct = pct/max(pct);

% Remove group delay due to filters so that impulse response is
% centered at zero.
tct = (0:length(pct)-1)*dt;
tct = tct - ceil(tx.gdac_delay + tx.hl_delay*ofdm.fsct + rx.gadc_delay)*dt;

n0 = ceil(tx.gdac_delay + tx.hl_delay*ofdm.fsct + rx.gadc_delay) + 1; % new zero after remove group delay
% Sampling at the chip rate
tn = tct([n0:-sim.Mct:1, n0+sim.Mct:sim.Mct:end]);
pn = pct([n0:-sim.Mct:1, n0+sim.Mct:sim.Mct:end]);

figure
subplot(211)
plot(tct, pct)
hold on
stem(tn, pn, 'fill')
plot(-Nneg_os/ofdm.fs*[1 1], [-1 1], 'k')
plot(Npos_os/ofdm.fs*[1 1], [-1 1], 'k')
text(-(Nneg_os)/ofdm.fs, 0.7, sprintf('N_{neg} = %d', Nneg_os))
text((Npos_os+0.1)/ofdm.fs, 0.7, sprintf('N_{pos} = %d', Npos_os))
axis([1.5*Npre_os*[-1 1]*1/ofdm.fs -0.5 1])
xlabel('t (s)')
ylabel('p(t)')
title('Impulse response of the channel (DAC * Laser * ADC)')

Npoints = 200;      % Number of points to plot
indplot = 1:ceil(length(f)/Npoints):length(f); % indexes of points to plot

% CD frequency response
theta = -1/2*D2beta2(fiber.D(tx.lamb), tx.lamb).*(2*pi*f(indplot)).^2*fiber.L; % theta = -1/2*beta2*w.^2*L
Hcdplot = cos(theta) - tx.alpha*sin(theta);        % fiber frequency response

subplot(212)
plot(f(indplot)/1e9, 20*log10(abs(Gdac(indplot))),...
f(indplot)/1e9, 20*log10(abs(Hl(indplot))),...
f(indplot)/1e9, 20*log10(abs(Gadc(indplot))),...
f(indplot)/1e9, 20*log10(abs(Hcdplot)),...
f(indplot)/1e9, 20*log10(abs(Gdac(indplot).*Hl(indplot).*Hcdplot.*Gadc(indplot))))
xlabel('f (GHz)')
ylabel('20log_{10}(H(f))')
axis([0 ofdm.fs/2e9 -50 5])
legend('DAC', 'Laser', 'ADC', 'Fiber', 'Total', 'Location', 'SouthWest')  

% 
% %% Transmitted signal
% Nsymb_plot = 50;        % Number of symbols to plot
% 
% indplot = 1:sim.Mct:sim.Mct*Nsymb_plot*(Nc+Npre_os); % Downsample before plotting
% tplot = sim.t(indplot);
% 
% figure
% subplot(221)
% plot(tplot, xtc(indplot))
% xlabel('t')
% ylabel('x_+(t)')
% axis([0 tplot(end) 0 1.2*max(xtc)])
% 
% subplot(222)
% plot(tplot, Pt(indplot)*1e3, ...
%     tplot([1 end]), Ptxest*[1 1]*1e3, 'k',...
%     tplot([1 end]), mean(Pt)*[1 1]*1e3, 'k')
% xlabel('t')
% ylabel('P(t) (mW)')
% axis([0 tplot(end) 0 1.2*max(Pt)*1e3])
% 
% subplot(223)
% spectrum_periodogram(AmpNorm(xtc(indplot) - mean(xtc)), Nc+Npre_os, ofdm.fs);
% 
% subplot(224)
% spectrum_periodogram(AmpNorm(Pt(indplot)-mean(Pt)), Nc+Npre_os, ofdm.fs);

%% Received signal vs transmitted
Npoints = min(2^10, length(dataRXm));         % number of points to show in the constellation

figure
subplot(221)
plot(dataRXm(1, 1:Npoints), '.')
hold on
plot(qammod(0:CS(1)-1, CS(1)), 'xk')
axis(sqrt(max(CS(1)))*[-1 1 -1 1])
axis square
title('1st subcarrier');

subplot(222)
plot(dataRXm(end, 1:Npoints), '.')
hold on
plot(qammod(0:CS(end)-1, CS(end)), 'xk')
axis(sqrt(max(CS(end)))*[-1 1 -1 1])
axis square
title('Last subcarrier');

subplot(223)
stem(1:Nu/2, log2(CS), 'filled')
xlabel('Subcarrier')
ylabel('log_2(CS)')

subplot(224)
stem(1:Nu/2, SNRn, 'fill')
xlabel('Subcarrier')
ylabel('SNR (dB)');

%% Quantization
if sim.quantiz
    tplot = sim.t(1:sim.Mct:100*sim.Mct*(Npre_os+Nc));
    figure
    subplot(321)
    plot(tplot, xncp(1:100*(Npre_os+Nc)), 'b')
    hold on
    plot(tplot([1 end]), yqtx([1 1]), 'k')
    plot(tplot([1 end]), yqtx([end end]), 'k')
    axis([tplot(1) tplot(end) 1.2*yqtx([1 end])])
    title('Input signal of quantiz (TX)')

    subplot(322)
    plot(tplot, yncp(1:100*(Npre_os+Nc)), 'b')
    hold on
    plot(tplot([1 end]), yqrx([1 1]), 'k')
    plot(tplot([1 end]), yqrx([end end]), 'k')
    axis([tplot(1) tplot(end) 1.2*yqrx([1 end])])
    title('Input signal of quantiz (RX)')

    subplot(323)
    hist(xncp, 50)
    hold on
    axis manual
    plot(yqtx([1 1]), [0 1e6], 'k')
    plot(yqtx([end end]), [0 1e6], 'k')
    title('Distribution of input signal (TX)')

    subplot(324)
    hist(yncp, 50)
    hold on
    axis manual
    plot(yqrx([1 1]), [0 1e6], 'k')
    plot(yqrx([end end]), [0 1e6], 'k')
    title('Distribution of input signal (RX)')

    subplot(325)
    hist(QuantErrTX, 2^sim.ENOB);
    title('Quantization noise (TX)')

    subplot(326)
    hist(QuantErrRX, 2^sim.ENOB);       
    title('Quantization noise (RX)')
end

% %% Clipping
% varDk = clip_noise_var(xncpc - K*xncp, ofdm, sigtx^2);

%% Number of errors
figure
subplot(121)
stem(1:Nu/2, numerr, 'fill')
xlabel('Subcarrier')
ylabel('Number of errors');
subplot(122)
[~, kk] = max(numerr);
plot(real(dataRXm(kk,:)), imag(dataRXm(kk,:)), '.')
hold on
plot(qammod(0:CS(kk)-1, CS(kk)), 'xk')
axis square
axis(sqrt(max(CS(kk)))*[-1 1 -1 1])
title(sprintf('Subcarrier #%d (worst subcarrier)', kk))
