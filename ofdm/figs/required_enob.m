%% Required ENOB to transmitt 16-QAM and 64-QAM for ACO-OFDM and DC-OFDM
clear, clc, close all

Nc = 64;
Nu = 52;
ros = Nc/Nu;

r = linspace(3, 5);

%% Anonymous function
% Require ENOB for DC-OFDM (SNR must be in linear units)
ENOBdcofdm = @(r, SNRn) 0.5*log2(2*r.^2/(3*ros).*SNRn);

% Require ENOB for ACO-OFDM (SNR must be in linear units)
% ENOBacoofdm = @(r, SNRn) 0.5*log2((r.^2 + 2*(1/sqrt(2*pi) + r).^2)/(6*ros).*SNRn);
ENOBacoofdm = @(r, SNRn) 0.5*log2((r.^2 + 2*(1/sqrt(2*pi) + r).^2)/(12*ros).*SNRn);


%% SNRn
M = 16;
SNRdB16 = fzero(@(x)berawgn(x - 10*log10(log2(M)),'qam',M) - 1.8e-4, 20);
SNR16 = 10^(SNRdB16/10);

M = 64;
SNRdB64 = fzero(@(x)berawgn(x - 10*log10(log2(M)),'qam',M)-1.8e-4, 20);
SNR64 = 10^(SNRdB64/10);


%% DC-OFDM
dcofdm_ENOB16 = ENOBdcofdm(r, SNR16);

dcofdm_ENOB64 = ENOBdcofdm(r, SNR64);

%% ACO-OFDM
acoofdm_ENOB16 = ENOBacoofdm(r, SNR16);

acoofdm_ENOB64 = ENOBacoofdm(r, SNR64);


%% Plot
figure1 = figure('Color',[1 1 1]);
hold on
plot(r, dcofdm_ENOB16, '--k', 'LineWidth', 1.5)
plot(r, acoofdm_ENOB16, '-k', 'LineWidth', 1.5)

plot(r, dcofdm_ENOB64, '--k', 'LineWidth', 1.5)
plot(r, acoofdm_ENOB64, '-k', 'LineWidth', 1.5)

legend('DC-OFDM', 'ACO-OFDM', 'Location', 'SouthEast');

% Create xlabel
xlabel('Clipping ratio','FontSize',14);

% Create ylabel
ylabel('Required ENOB','FontSize',14);

set(gca, 'FontSize', 14)
set(gca, 'xtick', 2:0.25:5)
axis([r(1) r(end) 3 6])
grid on
box on

% Create textbox
annotation(figure1,'textbox',...
    [0.75 0.44 0.14 0.06],...
    'String',{'16-QAM'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.75 0.85 0.14 0.06],...
    'String',{'64-QAM'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');


saveas(gca, 'required_enob', 'emf')

% convert gca to latex
% matlab2tikz files must be in current directory
matlab2tikz('tikz\required_enob.tikz')
