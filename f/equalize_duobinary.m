function [yd, eq] = equalize_duobinary(eq, yt, mpamdb, sim, verbose)      
%% Adaptive time-domain fractionally-spaced equalizer for duobinary M-PAM
%% This is only used if all duobinary levels are received. In the case that the decoding was already performed a typical equalizer can be used

% Check if number of taps was defined and is even. If even, make it odd
if isfield(eq, 'Ntaps') && mod(eq.Ntaps, 2) == 0
    eq.Ntaps = eq.Ntaps + 1;
end

transpose = false;
if size(yt, 1) < size(yt, 2)
    transpose = true;
    yt = yt.';
end

%% Filters are defined as follows
%      _____           _____
%     |     |    /    |     |
% --->| Hrx |---/ --->| Hff |--->
%     |_____|         |_____|
Ntaps = eq.Ntaps;
ros = eq.ros;
Ntrain = eq.Ntrain;
mu = eq.mu;
trainingSymbols = eq.trainSeq;

% Number of filters
[Nsamp, Nfilters]  = rat(ros);

% Initialize filter coefficients
W = zeros(Ntaps, Nfilters);
W(ceil(Ntaps/2), :) = 1;
Wdc = zeros(1, Nfilters); % adjustable DC bias

e = zeros(1, ceil(length(yt)/ros));
filt = 1;
window = (1:Ntaps) - ceil(Ntaps/2); % {-(Ntaps-1)/2,...,0,... (Ntaps-1)/2}.
kstart = ceil(Ntaps/(2*Nsamp))*Nsamp+1; % loop through detected symbols
k = kstart;
yd = zeros(1, sim.Nsymb); % output
for n = ((kstart-1)/ros+1):sim.Nsymb-(ceil(Ntaps/2)+sim.Ndiscard)
    % k indexes samples and n indexes symbols
    z = yt(k + window); % input of filter

    yd(n) = W(:, filt).'*z + Wdc(filt);

    if n < Ntrain % Training
        e(n) = yd(n) - trainingSymbols(n);
    else % decision directed
        e(n) = yd(n) - mpamdb.mod(mpamdb.demod(abs(yd(n))));
    end

    W(:, filt) = W(:, filt) - 2*mu*e(n)*z;
    Wdc(filt) = Wdc(filt) - 2*mu*e(n)*1;

    filt = mod(filt, Nfilters) + 1;

    % ensures that the center of the window of samples remains 
    % close to the nth symbol
    if abs((k+1)/ros - n-1) > 0.5
        k = k + 2;
    else
        k = k + 1;
    end
end

eq.h = W; % [Ntaps x Nfilters]
eq.Wdc = Wdc; % ajustable dc bias
% Frequency response is calculated for only one of the filters
eq.Hff = @(f) freqz(eq.h(:, 1), 1, 2*pi*f).*exp(1j*2*pi*f*grpdelay(eq.h(:, 1), 1, 1)); % removed group delay
eq.MSE = abs(e).^2;

if exist('verbose', 'var') && verbose       
    figure(400), clf
    subplot(221), hold on, box on
    plot(abs(e).^2)
    a = axis;
    plot(eq.Ndiscard(1)*[1 1], a(3:4), ':k')
    plot((sim.Nsymb-eq.Ndiscard(2))*[1 1], a(3:4), ':k')
    legend('MSE', 'Valid window')
    xlabel('Iteration')
    ylabel('MSE')
    title('Equalizer MSE')


    fnorm = linspace(-0.5, 0.5);
    for k = 1:size(W, 2)
        subplot(222), hold on, box on
        stem(-(eq.Ntaps-1)/2:(eq.Ntaps-1)/2, abs(W(:, k)))
        plot([-(eq.Ntaps-1)/2 (eq.Ntaps-1)/2], Wdc(k)*[1 1])

        Hw = freqz(W(:, k), 1, 2*pi*fnorm);
        subplot(223), hold on, box on
        plot(mpamdb.Rs*eq.ros*fnorm/1e9, abs(Hw).^2)

        subplot(224), hold on, box on
        plot(mpamdb.Rs*eq.ros*fnorm/1e9, unwrap(angle(Hw)))
    end
    subplot(222)   
    xlabel('n')
    ylabel('|W(n)|')
    title('Equalizer taps')
    legend('Equalizer taps', 'DC bias')

    subplot(223)
    a = axis;
    axis([-mpamdb.Rs*eq.ros/2e9 mpamdb.Rs*eq.ros/2e9 a(3:4)]);
    xlabel('Frequency (GHz)')
    ylabel('|H_W(f)|^2')

    subplot(224)
    a = axis;
    axis([-mpamdb.Rs*eq.ros/2e9 mpamdb.Rs*eq.ros/2e9 a(3:4)]);
    xlabel('Frequency (GHz)')
    ylabel('arg(H_W(f))')
    drawnow
end

yd = abs(yd);
if transpose
    yd = yd.';
end