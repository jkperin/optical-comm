function [yd, eq] = nonlinear_equalizer(eq, yt, mpam, sim, verbose)
%% Design equalizer and equalize yt
% If yt is empty only designs equalizer
% Inputs: 
% - eq: equalizer struct with equalizer parameters
%       - ros: oversampling ratio at which equalizer works. Typically the
%       same as receiver DSP ros.
%       - Ntaps: number of equalizer taps
%       - Hrx: receiver filter
%       - NSR: noise-to-signal ratio. If provided equalize.m designs MMSE
%       linear equalizer, otherwise zero-forcing liner equalizer is
%       designed
%       - Ntrain (if adaptive): training sequence length
%       - mu (if adaptive): adaptation rate
%       - MSE (if adaptive): mean square error
% - yt: input signal
% - Hch (required if fixed equalizer): channel frequency response including pulse shape
% - mpam: mpam class
% - sim: simulation parameters struct
% - verbose (optional, default = false): whether to plot results
% Outputs:
% - yd: equalized signal sampled at symbol rate
% - eq: equalizer struct with equalizer parameters
% These are the fields added by equalize.m
%       - Hff: function handle for the frequency response of equalizer
%       - h: taps of equalizer (Ntaps x Nfilters)

% Equalizer types supported:
% - None: No equalizer. Just sample at symbol rate.
% - Fixed or Adaptive TD-SR-LE: time-domain symbol-rate linear equalizer.
% This assumes that there is a matched filter before symbol-rate sampling.
% - Adaptive TD-FS-LE: time-domain fractionally-spaced linear equalizer.
% This assumes an anti-aliasing filter before sampling.

% Check if number of taps was defined and is even. If even, make it odd
if isfield(eq, 'Ntaps') && mod(eq.Ntaps, 2) == 0
    eq.Ntaps = eq.Ntaps + 1;
end

% Makes sure everything is in right dimensions
if size(sim.f, 1) == 1
    sim.f = sim.f.';
end

ytsize = size(yt);
if ytsize(1) == 1
    yt = yt.';
end

%% Adaptive nonlinear time-domain fractionally-spaced equalizer
Ntaps = eq.Ntaps;
ros = eq.ros;
Ntrain = eq.Ntrain;
mu = eq.mu;
trainingSymbols = mpam.mod(eq.trainSeq);

idx = nchoosek(1:Ntaps, 2);
Len = 2*Ntaps + size(idx, 1);
% Len = Ntaps;

% Number of filters
[Nsamp, Nfilters]  = rat(ros);

% Initialize filter coefficients
W = zeros(Len, Nfilters);
W(ceil(Len/2), :) = 1;
Wdc = zeros(1, Nfilters); % adjustable DC bias

e = zeros(1, ceil(length(yt)/ros));
filt = 1;
window = (1:Ntaps) - ceil(Ntaps/2); % {-(Ntaps-1)/2,...,0,... (Ntaps-1)/2}.
kstart = ceil(Ntaps/(2*Nsamp))*Nsamp+1; % loop through detected symbols
k = kstart;
yd = zeros(1, sim.Nsymb); % output
for n = ((kstart-1)/ros+1):sim.Nsymb-(ceil(Ntaps/2)+sim.Ndiscard)
    % k indexes samples and n indexes symbols
    z = yt(k + window); % input signal
    
    % Expand with nonlinear terms
    u = [z; z.^2; z(idx(:, 1)).*z(idx(:, 2))];
%     u = z;

    yd(n) = W(:, filt).'*u + Wdc(filt);

    if n < Ntrain % Training
        e(n) = yd(n) - trainingSymbols(n);
    else % decision directed
        e(n) = yd(n) - mpam.mod(mpam.demod(yd(n)));
    end

    W(:, filt) = W(:, filt) - 2*mu*e(n)*u;
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
    subplot(211), hold on, box on
    plot(abs(e).^2)
    a = axis;
    plot(eq.Ndiscard(1)*[1 1], a(3:4), ':k')
    plot((sim.Nsymb-eq.Ndiscard(2))*[1 1], a(3:4), ':k')
    legend('MSE', 'Valid window')
    xlabel('Iteration')
    ylabel('MSE')
    title('Equalizer MSE')

    for k = 1:size(W, 2)
        subplot(212), hold on, box on
        stem(-(Len-1)/2:(Len-1)/2, abs(W(:, k)))
        plot([-(Len-1)/2 (Len-1)/2], Wdc(k)*[1 1])
    end 
    xlabel('n')
    ylabel('|W(n)|')
    title('Equalizer taps')
    legend('Equalizer taps', 'DC bias')
end


% returns yd with right dimensions
if all(ytsize ~= size(yd))
    yd = yd.';
end