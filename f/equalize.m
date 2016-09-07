function [yd, eq] = equalize(eq, yt, Hch, mpam, sim, verbose)
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
    Hch = Hch.';
end

ytsize = size(yt);
if ytsize(1) == 1
    yt = yt.';
end

% if yt is empty i.e., only design filter, then change equalization type
% from adaptive to fixed
if isempty(yt)
    eq.type = strrep(eq.type, 'Adaptive', 'Fixed');
end

%% Filters are defined as follows
%      _____           _____
%     |     |    /    |     |
% --->| Hrx |---/ --->| Hff |--->
%     |_____|         |_____|

switch lower(eq.type)
    case 'none'
        %% No equalization
        eq.Hff = @ (f) ones(size(f)); % equalizer
        eq.h = 1;
        assert(floor(eq.ros) == ceil(eq.ros), 'equalize: when eq.type = none, eq.ros must be integer');
        yd = yt(1:eq.ros:end);
               
    case 'adaptive td-le'
        %% Adaptive Time-domain fractionally-spaced equalizer
        Ntaps = eq.Ntaps;
        ros = eq.ros;
        Ntrain = eq.Ntrain;
        mu = eq.mu;
        trainingSymbols = mpam.mod(eq.trainSeq);
        
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
                e(n) = yd(n) - mpam.mod(mpam.demod(yd(n)));
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
                plot(mpam.Rs*eq.ros*fnorm/1e9, abs(Hw).^2)

                subplot(224), hold on, box on
                plot(mpam.Rs*eq.ros*fnorm/1e9, unwrap(angle(Hw)))
            end
            subplot(222)   
            xlabel('n')
            ylabel('|W(n)|')
            title('Equalizer taps')
            legend('Equalizer taps', 'DC bias')
            
            subplot(223)
            a = axis;
            axis([-mpam.Rs*eq.ros/2e9 mpam.Rs*eq.ros/2e9 a(3:4)]);
            xlabel('Frequency (GHz)')
            ylabel('|H_W(f)|^2')

            subplot(224)
            a = axis;
            axis([-mpam.Rs*eq.ros/2e9 mpam.Rs*eq.ros/2e9 a(3:4)]);
            xlabel('Frequency (GHz)')
            ylabel('arg(H_W(f))')
            drawnow
        end
        
    case 'fixed td-sr-le' 
        if not(isfield(eq, 'ros'))
            eq.ros = 1;
        elseif eq.ros ~= 1
            warning('equalize/eq.ros = %.2f, but symbol-rate equalization is set', eq.ros)
        end
        Ntaps = eq.Ntaps;
        
        Hch = Hch/sqrt(interp1(sim.f, abs(Hch).^2, 0));   % normalize to unit gain at DC 
        Hmatched = conj(Hch); % matched filterd matched to the received pulse shape

        %% MMSE Time-domain symbol-rate equalizer
        n = -floor(Ntaps/2)*sim.Mct:sim.Mct*floor(Ntaps/2);

        hmatched2 = real(ifft(ifftshift(Hmatched)));
        hmatched = hmatched2(1:sim.Mct*floor(Ntaps/2));
        hmatched = [hmatched2(end-sim.Mct*floor(Ntaps/2):end); hmatched];

        x = conv(hmatched, hmatched(end:-1:1), 'same');
        xd = x(mod(abs(n), sim.Mct) == 0); % symbol-rate sample received pulse

        xd = xd/abs(sum(xd)); % normalize to unit gain at DC

        X = toeplitz([xd; zeros(Ntaps-1, 1)], [xd(1) zeros(1, Ntaps-1)]);

        % Get correct rows from Toeplitz matrix
        X = X(ceil(size(X, 1)/2)+(-floor(Ntaps/2):floor(Ntaps/2)), :);

        % ZF condition
        e = zeros(Ntaps, 1); 
        e((Ntaps+1)/2) = 1;

        % NSR = noise signal ratio
        if isfield(eq, 'NSR')
            %% MMSE linear equalizer
            W = X*(((X' + eq.NSR*eye(Ntaps))*X)\e);
        else
            %% ZF linear equalizer
            W = X*((X'*X)\e);
        end                
        % Note: if NSR = 0, or if NSR is not specified, then
        % zero-forcing equalizer is designed 

        % Filter
        if isempty(yt) % only design
            yd = []; 
        else % filter
            yd = filter(W, 1, yt);
            yd = circshift(yd, [-(Ntaps-1)/2 0]); % remove delay due to equalizer
        end  

        eq.hmatched = hmatched;

        % Aux
        eq.h = W;
        eq.Hff = @(f) freqz(eq.h, 1, 2*pi*f).*exp(1j*2*pi*f*grpdelay(eq.h, 1, 1)); % removed group delay
        eq.Hrx = Hmatched.'; 
        
        if exist('verbose', 'var') && verbose       
            figure(400)
            subplot(211), box on
            stem(-(eq.Ntaps-1)/2:(eq.Ntaps-1)/2, abs(W))
            xlabel('n')
            ylabel('|W(n)|')
            title('Equalizer taps')
            subplot(212), box on
            fnorm = linspace(-0.5, 0.5);
            Hw = freqz(W, 1, 2*pi*fnorm);
            plot(fnorm*mpam.Rs*eq.ros/1e9, abs(Hw).^2)
            a = axis;
            axis([-mpam.Rs*eq.ros/2e9 mpam.Rs*eq.ros/2e9 a(3:4)]);
            xlabel('Frequency (GHz)')
            ylabel('|H_W(f)|^2')
            title('Equalizer frequency response')
            drawnow
        end
    otherwise
        error('equalize: Equalization type *%s* not implemented yet!', eq.type)
end

% returns yd with right dimensions
if all(ytsize ~= size(yd))
    yd = yd.';
end