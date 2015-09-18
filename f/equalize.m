%% Design equalizer and equalize yt
function [yd, eq] = equalize(eq, yt, Hch, mpam, rx, sim)
% Types supported:
% - None: No equalizer. Filter using antialiasing filter (rx.elefilt) and
% sample at symbol rate.
% - Analog: Invert channel response and apply matched-filter matched to the
% pulse shape and then do symbol-rate sampling
% - Fixed or Adaptive TD-SR-LE: time-domain symbol-rate linear equalizer.
% This assumes that there is a matched filter before symbol-rate sampling.
% - Fixed or Adaptive TD-FS-LE: time-domain fractionally-spaced linear equalizer.
% This assumes an anti-aliasing filter before sampling.

% Check if number of taps was defined and is even. If even, make it odd
if isfield(eq, 'Ntaps') && mod(eq.Ntaps, 2) == 0
    eq.Ntaps = eq.Ntaps + 1;
end

% If channel response isn't define don't do equalization
if isempty(Hch)
    eq.type = 'None';
else % normalize channel response to have unit gain at DC
    Gtx = design_filter('matched', mpam.pshape, 1/sim.Mct); % transmitted pulse shape
    
    Gtx = Gtx.H(sim.f/sim.fs);

    Hch = Hch/interp1(sim.f, Hch, 0);    
end

% if yt is empty i.e., only design filter, then change equalization type
% from adaptive to fixed
if isempty(yt)
    eq.type = strrep(eq.type, 'Adaptive', 'Fixed');
end

% Inserts delay of half of sample if oversampling of continuous time is even
if mod(sim.Mct, 2) == 0 
    Delay = exp(-1j*pi*sim.f/sim.fs);
else 
    Delay = 1;
end

switch eq.type
    case 'None'
        %% No equalization
        eq.H = ones(size(sim.f));
        eq.f = sim.f;
        eq.Kne = 1;
               
        if ~isempty(yt)
            % Receiver filter
            yt = real(ifft(fft(yt).*ifftshift(Delay.*rx.elefilt.H(sim.f/sim.fs))));

            yd = yt(floor(sim.Mct/2)+1:sim.Mct:end); % +1 because indexing starts at 1
        else
            yd = [];
        end
    case 'Analog'
        %% Analog Equalization
        % -> Channel inversion filter -> matched filter matched to
        % transmitted pulse.
        % Channel inversion filter is defined as follows
        % 1/Channel for |f| < mpam.Rs
        % 0 for |f| > mpam.Rs   
        
        Heq = zeros(size(sim.f));
        Heq(abs(sim.f)<mpam.Rs) = 1./(Hch(abs(sim.f)<mpam.Rs));
        
        % Channel inversion + Matched filter
        Grx = Heq.*Delay.*conj(Gtx.*Heq);

        % Apply equalization filter
        if isempty(yt) % function was called just to calculate equalizer
            yd = [];
        else
            % Channel inversion + Matche filter
            yteq = real(ifft(fft(yt).*ifftshift(Grx))); 
            
            % Sample
            yd = yteq(floor(sim.Mct/2)+1:sim.Mct:end); % +1 because indexing starts at 1
        end
        
        eq.H = Heq;
        eq.f = sim.f;
        
        eq.Kne = trapz(sim.f, abs(Grx).^2)/(mpam.Rs);
        
    case 'Fixed TD-FS-LE'
        % include transmitted pulse shape and receiver antialiasing filter
        Hch = Gtx.*Hch.*rx.elefilt.H(sim.f/sim.fs); 
        
        Hmatched = Delay.*conj(Hch);

        %% MMSE Time-domain symbol-rate equalizer
        n = -floor(eq.Ntaps/2)*sim.Mct:sim.Mct*floor(eq.Ntaps/2);
        nd = [fliplr(-sim.Mct/eq.ros:-sim.Mct/eq.ros:n(1)) 0:sim.Mct/eq.ros:n(end)];

        % Design matched filter           
        hmatched2 = real(ifft(ifftshift(Hmatched)));
        hmatched = hmatched2(1:sim.Mct*floor(eq.Ntaps/2));
        hmatched = [hmatched2(end-sim.Mct*floor(eq.Ntaps/2):end); hmatched];

        % p = matched filter, x = pulse shape after matched filter
        x = conv(hmatched, hmatched(end:-1:1), 'same');
        xd = interp1(n, x, nd, 'spline'); % oversampled impulse response
        pd = interp1(n, hmatched, nd, 'spline'); % oversampled impulse response
        pd = pd(ceil(length(pd)/2)+(-floor(eq.Ntaps/2):floor(eq.Ntaps/2)));
        xd = xd/abs(sum(xd)); % normalize to unit gain at DC
        pd = pd/abs(sum(pd)); % normalize to unit gain at DC

        % Design equalizer
        X = toeplitz([xd.'; zeros(eq.Ntaps-1, 1)], [xd(1) zeros(1, eq.Ntaps-1)]);

        % Get correct rows from Toeplitz matrix
        X = X(ceil(size(X, 1)/2)+(-floor(eq.Ntaps/2):floor(eq.Ntaps/2)), 1:eq.ros:end);

        e = zeros((eq.Ntaps+1)/2, 1); 
        e((length(e)+1)/2) = 1;

        % NSR = noise signal ratio
        if isfield(eq, 'NSR')
            W = X*(((X' + eq.NSR*eye(eq.Ntaps))*X)\e);
        else
            W = X*((X'*X)\e);
        end                
        % Note: if NSR = 0, or if NSR is not specified, then
        % zero-forcing equalizer is designed    

        % Convolve with matched filter
        W = conv(W, pd, 'same');
        W = W/abs(sum(W)); % normalize to unit gain at DC

        % Filter using
        if isempty(yt) % only design
            yd = []; 
        else
            % Antialiasing filter
            ytaa = real(ifft(fft(yt).*ifftshift(rx.elefilt.H(sim.f/sim.fs))));
            
            if mod(sim.Mct/eq.ros, 2) == 0
                yk = ytaa(1:sim.Mct/eq.ros:end);
%                 tk = sim.t(1:sim.Mct/eq.ros:end);
            else % if sim.Mct is not multiple of ros, then interpolate
                yk = interp1(1:length(ytaa), ytaa, 1:sim.Mct/eq.ros:length(ytaa), 'spline');
%                 tk = interp1(1:length(ytaa), sim.t, 1:sim.Mct/eq.ros:length(ytaa));

    %             yk = resample(ytaa, sim.ros, sim.Mct);
    %             tk = resample(sim.t, sim.ros, sim.Mct);            
            end
            % Filter
            yk = filter(W, 1, yk); 
            yk = circshift(yk, [0 -(eq.Ntaps-1)/2]); % remove delay of FIR filter
            % Note: !! In this case, W is not necessarily linear phase, so
            % the delay (eq.Ntaps-1)/2 is not necessarily exact
            yd = yk(2:eq.ros:end).';
        end  

        % Aux
        eq.num = W;
        eq.den = 1;
        [eq.H, w] = freqz(eq.num, eq.den);
        eq.f = w/(2*pi);
        eq.Kne = 2*eq.ros*trapz(eq.f, abs(eq.H));
        % Note: eq.f is one-sided (x2) and the sampling rate here is rosxRs
        
    case 'Adaptive TD-FS-LE'
        %% Adaptive Time-domain fractionally-spaced equalizer
        if mod(rx.eq.Ntaps, 2) == 0
            eq.Ntaps = eq.Ntaps + 1;
        end
        Ntaps = eq.Ntaps;
        ros = eq.ros;
        Ntrain = eq.Ntrain;
        mu = eq.mu;
        b = mpam.mod(rx.eq.TrainSeq, 1);
        
        % Antialiasing filter
        ytaa = real(ifft(fft(yt).*ifftshift(Delay.*rx.elefilt.H(sim.f/sim.fs))));
        
        if mod(sim.Mct/ros, 2) == 0
            yk = ytaa(1:sim.Mct/ros:end);
            tk = sim.t(1:sim.Mct/ros:end);
        else % if sim.Mct is not multiple of ros, then interpolate
            yk = interp1(1:length(ytaa), ytaa, 1:sim.Mct/ros:length(ytaa));
            tk = interp1(1:length(ytaa), sim.t, 1:sim.Mct/ros:length(ytaa));
            
%             yk = resample(ytaa, sim.ros, sim.Mct);
%             tk = resample(sim.t, sim.ros, sim.Mct);            
        end
        
        if size(yk, 1) < size(yk, 2)
            yk = yk.';
        end
                
        W = zeros(Ntaps, 1); % Filter taps
        W((Ntaps+1)/2) = 1;
        y = zeros(size(yk));
        n = 1;
        e = zeros(1, length(yk)/ros);
        for k = max(Ntaps, sim.Ndiscard*ros):length(yk)
            z = yk(k-Ntaps+1:k);

            y(k) = sum(W.*z);

            if mod(k-(Ntaps+1)/2, ros) == 0
                if n < Ntrain % Training
                    e(k/ros) = y(k) - b((k-(Ntaps+1)/2)/ros);
                else
                    e(k/ros) = y(k) - mpam.mod(mpam.demod(y(k)), 1);
                end
                n = n + 1;

                W = W - 2*mu*e(k/ros)*z;
            end
        end
        
        % remove delay inserted by the transversal FIR filter
        y = circshift(y, [-(Ntaps+1)/2 0]);

        yd = y(ros:ros:end);
        
        if sim.verbose
            figure, plot(e)
            xlabel('Iteration')
            ylabel('Error')
            
            figure, hold on
            plot(sim.t, yt)
            plot(tk, yk, 'o')
            plot(tk(ros:ros:end), yd, '*')
            
            [H, w] = freqz(W(end:-1:1), 1);
            figure, hold on
            plot(2*mpam.Rs/1e9*w/(2*pi), abs(H).^2)
        end
        
        eq.num = W(end:-1:1);
        eq.den = 1;
        [eq.H, w] = freqz(eq.num, eq.den);
        eq.f = w/(2*pi);
        eq.Kne = 2*eq.ros*trapz(eq.f, abs(eq.H).^2);
        % Note: eq.f is one-sided (x2) and the sampling rate here is rosxRs
        
    case 'Adaptive TD-SR-LE'
        Ntaps = eq.Ntaps;
        Ntrain = eq.Ntrain;
        mu = eq.mu;
        b = mpam.mod(rx.eq.TrainSeq, 1);
        
        % Received Pulse Spectrum
        Hmatched = Delay.*conj(Gtx.*Hch);

        yt = real(ifft(fft(yt).*ifftshift(Hmatched)));

        yk = yt(floor(sim.Mct/2)+1:sim.Mct:end);
        
        W = zeros(Ntaps, 1);
        W((Ntaps+1)/2) = 1;
        yd = zeros(size(yk));
        e = zeros(1, length(yk));
        n = 1;
        for k = Ntaps:length(yk)
            z = yk(k-Ntaps+1:k);

            yd(k) = sum(W.*z);

            if n < Ntrain % Training
                e(k) = yd(k) - b(k-(Ntaps+1)/2);
            else
                e(k) = yd(k) - mpam.mod(mpam.demod(yd(k)), 1);
            end
            n = n + 1;

            W = W - 2*mu*e(k)*z;
        end

        % remove delay inserted by the transversal FIR filter
        yd = circshift(yd, [-(Ntaps+1)/2 0]);
        
%         figure, plot(yd, 'o')
%         figure, plot(e)
%         
        eq.num = W(end:-1:1);
        eq.den = 1;
        [eq.H, w] = freqz(eq.num, eq.den);
        eq.f = w/(2*pi);
        eq.Kne = 2*trapz(eq.f, abs(eq.H));

        case 'Fixed TD-SR-LE'            
            Hmatched = Delay.*conj(Gtx.*Hch);
            % Note: Fiber frequency response is real, thus its group delay
            % is zero.
                        
            %% MMSE Time-domain symbol-rate equalizer
            n = -floor(eq.Ntaps/2)*sim.Mct:sim.Mct*floor(eq.Ntaps/2);

            hmatched2 = real(ifft(ifftshift(Hmatched)));
            hmatched = hmatched2(1:sim.Mct*floor(eq.Ntaps/2));
            hmatched = [hmatched2(end-sim.Mct*floor(eq.Ntaps/2):end); hmatched];

            x = conv(hmatched, hmatched(end:-1:1), 'same');
            xd = x(mod(abs(n), sim.Mct) == 0); % symbol-rate sample received pulse

            xd = xd/abs(sum(xd)); % normalize to unit gain at DC
            
            X = toeplitz([xd; zeros(eq.Ntaps-1, 1)], [xd(1) zeros(1, eq.Ntaps-1)]);
            
            % Get correct rows from Toeplitz matrix
            X = X(ceil(size(X, 1)/2)+(-floor(eq.Ntaps/2):floor(eq.Ntaps/2)), :);
            
            % ZF condition
            e = zeros(eq.Ntaps, 1); 
            e((eq.Ntaps+1)/2) = 1;
            
            % NSR = noise signal ratio
            if isfield(eq, 'NSR')
                W = X*(((X' + eq.NSR*eye(eq.Ntaps))*X)\e);
            else
                W = X*((X'*X)\e);
            end                
            % Note: if NSR = 0, or if NSR is not specified, then
            % zero-forcing equalizer is designed 
            
            % Filter using
            if isempty(yt) % only design
                yd = []; 
            else
                % Matched filter
                yt = real(ifft(fft(yt).*ifftshift(Hmatched)));
                % Sample at symbol rate
                yk = yt(floor(sim.Mct/2)+1:sim.Mct:end);
                % Filter
                yd = filter(W, 1, yk);
                yd = circshift(yd, [-(eq.Ntaps-1)/2 0]); % remove delay due to equalizer
            end  

            % Aux
            eq.num = W;
            eq.den = 1;
            [eq.H, w] = freqz(eq.num, eq.den);
            eq.f = w/(2*pi);
            eq.Kne = 2*trapz(eq.f, abs(eq.H));        
    otherwise
        error('Equalization option not implemented yet!')
end

if sim.verbose
    plot(sim.f/1e9, abs(Heq).^2)
    xlabel('Frequency (GHz)')
    ylabel('|H_{eq}(f)|^2')
    axis([0 2*mpam.Rs/1e9 0 1.2*max(abs(Heq).^2)])
end
    

