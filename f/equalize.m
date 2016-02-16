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

    Hch = Hch/interp1(sim.f, Hch, 0);   % normalize to unit gain at DC 
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


switch eq.type
    case 'None'
        %% No equalization
        eq.Kne = 1;
        eq.Hrx = rx.elefilt.H(sim.f/sim.fs); % receiver filter
        eq.Hff = @ (f) ones(size(f)); % equalizer
               
        if ~isempty(yt)
            % Receiver filter
            yt = real(ifft(fft(yt).*ifftshift(eq.Hrx)));

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
        Grx = Heq.*conj(Gtx.*Heq);

        % Apply equalization filter
        if isempty(yt) % function was called just to calculate equalizer
            yd = [];
        else
            % Channel inversion + Matche filter
            yteq = real(ifft(fft(yt).*ifftshift(Grx))); 
            
            % Sample
            yd = yteq(floor(sim.Mct/2)+1:sim.Mct:end); % +1 because indexing starts at 1
        end
        
        eq.Hrx = Grx;
        eq.Hff = @(f) ones(size(f));
        
        eq.Kne = trapz(sim.f, abs(Grx).^2)/(mpam.Rs);
        
    case 'Fixed TD-LE'
        %% Fixed Time-domain fractionally-spaced equalizer
        ros = eq.ros;
        
        if mod(rx.eq.Ntaps, 2) == 0
            eq.Ntaps = eq.Ntaps + 1;
        end
        Ntaps = eq.Ntaps;
        
        % Number of filters
        [Nsamp, Nfilters]  = rat(ros);
       
        % Antialiasing / matched filter
        if ros ~= 1 % if oversampling use rx.elefilt as Antialiasing
            fprintf('> Antialiasing filter type: %s\n', rx.elefilt.type)
            Hrx = rx.elefilt.H(sim.f/sim.fs);       
        else % in case of symbol-rate sampling uses matched filter as receiver filter
            Hrx = conj(Gtx.*Hch);
            disp('> Antialiasing filter type: matched filter')
        end
        
        ytaa = real(ifft(fft(yt).*ifftshift(Hrx.*...
            exp(1j*2*pi*sim.f/sim.fs*(sim.Mct-1)/2)))); % time-shift signal
              % so that first sample corresponds to first symbol
        
        % Downsample
        if mod(sim.Mct/ros, 2) == 0
            yk = ytaa(1:sim.Mct/ros:end);
        else % if sim.Mct is not multiple of ros, then interpolate
            warning('equalize: sequence had to be interpolated because oversmapling ratio to simulation continuous time (Mct) is not a interger multiple of eq.ros')
            yk = interp1(1:length(ytaa), ytaa, 1:sim.Mct/ros:length(ytaa));   
        end    

        % Initialize filter coefficients
        W = zeros(Ntaps, Nfilters);
       
        filt = 1;
        window = (1:Ntaps) - ceil(Ntaps/2); % {-(Ntaps-1)/2,...,0,... (Ntaps-1)/2}.
        kstart = ceil(Ntaps/(2*Nsamp))*Nsamp+1; % loop through detected symbols
        k = kstart;
        yd = zeros(ceil(length(yk)/ros), 1); % output
        for n = ((kstart-1)/ros+1):length(yd)-(ceil(Ntaps/2)+sim.Ndiscard)
            % k indexes samples and n indexes symbols
            z = yk(k + window); % input of filter

            yd(n) = W(:, filt).'*z;

            filt = mod(filt, Nfilters) + 1;
            
            % ensures that the center of the window of samples remains 
            % close to the nth symbol
            if abs((k+1)/ros - n-1) > 0.5
                k = k + 2;
            else
                k = k + 1;
            end
        end
        
        Htot = Gtx.*Hch.*Hrx;

        %% MMSE Time-domain symbol-rate equalizer
        t = -floor(eq.Ntaps/2)/ros:1/sim.Mct:floor(eq.Ntaps/2)/ros;
        
        tmp = real(ifft(ifftshift(Htot)));
        htot = tmp(1:floor(length(t)/2));
        htot = [tmp(end-floor(length(t)/2):end); htot];
        
        % Initialize filter coefficients
        h = zeros(Ntaps, Nfilters);
        W = zeros(7, Nfilters);
        td = -floor(eq.Ntaps/2)/ros:1/ros:floor(eq.Ntaps/2)/ros;
               
        for filt = 1:Nfilters
            h(:, filt) = interp1(t, htot, td - filt + 1).';
            
            h(:, filt) = h(:, filt)/abs(sum(h(:, filt)));
            
            X = toeplitz([h(:, filt); zeros(Ntaps-1, 1)], [h(1, filt) zeros(1, Ntaps-1)]);
            
            % Get correct rows from Toeplitz matrix
            X = X(Ntaps + (-floor(Ntaps/2):Nsamp:floor(Ntaps/2)), :);

            % No ISI condition
            e = zeros(Ntaps, 1); 
            e(ceil(Ntaps/2)) = 1;
            
            % Design filter
            if isfield(eq, 'SNR') % MMSE
                W(:, filt) = X*(((X' + 1/eq.SNR*eye(eq.Ntaps))*X)\e);
            else % ZF
                W(:, filt) = X*((X'*X)\e);
            end                
            % Note: if NSR = 0, or if NSR is not specified, then
            % zero-forcing equalizer is designed
        end

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

        eq.hmatched = hmatched;

        % Aux
        eq.num = W;
        eq.den = 1;
        eq.Hff = @(f) freqz(eq.num, eq.den, 2*pi*f).*exp(1j*2*pi*f*grpdelay(eq.num, eq.den, 1)); % removed group delay
        eq.Hrx = Hmatched;

        [eq.H, w] = freqz(eq.num, eq.den);
        eq.f = w/(2*pi);
        eq.Kne = 2*trapz(eq.f, abs(eq.H));     
        
    case 'Adaptive TD-LE'
        %% Adaptive Time-domain fractionally-spaced equalizer
        ros = eq.ros;
        Ntrain = eq.Ntrain;
        mu = eq.mu;
        b = mpam.mod(rx.eq.TrainSeq, 1);
        
        if mod(rx.eq.Ntaps, 2) == 0
            eq.Ntaps = eq.Ntaps + 1;
        end
        Ntaps = eq.Ntaps;
        
        % Number of filters
        [Nsamp, Nfilters]  = rat(ros);

        % Antialiasing / matched filter
        if ros ~= 1 % if oversampling use rx.elefilt as Antialiasing
            fprintf('> Antialiasing filter type: %s\n', rx.elefilt.type)
            ytaa = real(ifft(fft(yt).*ifftshift(rx.elefilt.H(sim.f/sim.fs).*...
                exp(1j*2*pi*sim.f/sim.fs*(sim.Mct-1)/2)))); % time-shift signal
            % so that first sample corresponds to first symbol
        else % in case of symbol-rate sampling uses matched filter as receiver filter
            Hmatched = conj(Gtx.*Hch);
            disp('> Antialiasing filter type: matched filter')
            ytaa = real(ifft(fft(yt).*ifftshift(Hmatched.*...
                exp(1j*2*pi*sim.f/sim.fs*(sim.Mct-1)/2)))); % time-shift signal  
        end
        
        % Downsample
        if mod(sim.Mct/ros, 2) == 0
            yk = ytaa(1:sim.Mct/ros:end);
            tk = sim.t(1:sim.Mct/ros:end);
        else % if sim.Mct is not multiple of ros, then interpolate
            warning('equalize: sequence had to be interpolated because oversmapling ratio to simulation continuous time (Mct) is not a interger multiple of eq.ros')
            yk = interp1(1:length(ytaa), ytaa, 1:sim.Mct/ros:length(ytaa));
            tk = interp1(1:length(ytaa), sim.t, 1:sim.Mct/ros:length(ytaa));     
        end
        
        if size(yk, 1) < size(yk, 2)
            yk = yk.';
        end   
        
        % Initialize filter coefficients
        W = zeros(Ntaps, Nfilters);
        W(ceil(Ntaps/2), :) = 1;
       
        e = zeros(1, ceil(length(yk)/ros));
        filt = 1;
        window = (1:Ntaps) - ceil(Ntaps/2); % {-(Ntaps-1)/2,...,0,... (Ntaps-1)/2}.
        kstart = ceil(Ntaps/(2*Nsamp))*Nsamp+1; % loop through detected symbols
        k = kstart;
        yd = zeros(size(b)); % output
        W_old = zeros(size(W));
        for n = ((kstart-1)/ros+1):length(b)-(ceil(Ntaps/2)+sim.Ndiscard)
            % k indexes samples and n indexes symbols
            z = yk(k + window); % input of filter

            yd(n) = W(:, filt).'*z;

            if n < Ntrain % Training
                e(n) = yd(n) - b(n);
            else
                e(n) = yd(n) - mpam.mod(mpam.demod(yd(n)), 1);
            end

            W(:, filt) = W(:, filt) - 2*mu*e(n)*z;
            
            filt = mod(filt, Nfilters) + 1;
            
            % ensures that the center of the window of samples remains 
            % close to the nth symbol
            if abs((k+1)/ros - n-1) > 0.5
                k = k + 2;
            else
                k = k + 1;
            end
        end
        
        figure(100), plot(abs(e).^2)
        xlabel('Iteration')
        ylabel('MSE')
        axis([1 length(e) 0 2*max(e(end-100:end))])
               
        if sim.verbose         
            figure, hold on
            plot(sim.t, yt)
            plot(tk, yk, 'o')
            plot(tk(1:ros:end), yd, '*')
            axis([tk(end-200) tk(end) -0.5 1.5])
            
            [H, w] = freqz(W(end:-1:1), 1);
            figure, hold on
            plot(w/(2*pi), abs(H).^2)
        end
               
        eq.num = W(end:-1:1);
        eq.den = 1;
        eq.Hff = @(f) freqz(eq.num, eq.den, 2*pi*f).*exp(1j*2*pi*f*grpdelay(eq.num, eq.den, 1)); % removed group delay
   
    case 'Fixed TD-SR-LE'            
        Hmatched = conj(Gtx.*Hch);
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
            yt = real(ifft(fft(yt).*ifftshift(Hmatched.*...
                exp(1j*2*pi*sim.f/sim.fs*(sim.Mct-1)/2)))); % time-shift signal
            % so that first sample corresponds to first symbol
            % Sample at symbol rate
            yk = yt(1:sim.Mct:end);
            % Filter
            yd = filter(W, 1, yk);
            yd = circshift(yd, [-(eq.Ntaps-1)/2 0]); % remove delay due to equalizer
        end  

        eq.hmatched = hmatched;

        % Aux
        eq.num = W;
        eq.den = 1;
        eq.Hff = @(f) freqz(eq.num, eq.den, 2*pi*f).*exp(1j*2*pi*f*grpdelay(eq.num, eq.den, 1)); % removed group delay
        eq.Hrx = Hmatched;

        [eq.H, w] = freqz(eq.num, eq.den);
        eq.f = w/(2*pi);
        eq.Kne = 2*trapz(eq.f, abs(eq.H));     
        
   case 'Fixed FD-FS-LE'   
        ros = eq.ros;       

        yin = yt;
        yt = real(ifft(fft(yt).*ifftshift(rx.elefilt.H(sim.f/sim.fs).*...
           exp(1j*2*pi*sim.f/sim.fs*(sim.Mct-1)/2))));
        xk = yt(1:sim.Mct/ros:end);
        
        df = 1/length(xk);
        f = (-0.5:df:0.5-df)';
        
        Hchin = Hch;
        Hch = interp1(sim.f/sim.fs, Hch, f/ros);
        
        Heq = 1./Hch;
        gd = -diff(unwrap(phase(Heq)))/df;
        gd = fftshift(gd);
        gd0 = gd(1);
                
        % Filter
        yk = real(ifft(fft(xk).*ifftshift(Heq.*exp(1j*2*pi*f*gd0))));
        
        % Downsample
        yd = downsample(yk, ros);

        % Aux
        1;
        eq.Hff = @(f) freqz(eq.num, eq.den, 2*pi*f); % must remove group delay
        eq.Hrx = rx.elefilt.H(sim.f/sim.fs);

        [eq.H, w] = freqz(eq.num, eq.den);
        eq.f = w/(2*pi);
        eq.Kne = 2*trapz(eq.f, abs(eq.H));     
    otherwise
        error('Equalization option not implemented yet!')
end

if sim.verbose
    plot(sim.f/1e9, abs(eq.Hrx).^2)
    xlabel('Frequency (GHz)')
    ylabel('|H_{eq}(f)|^2')
    axis([0 2*mpam.Rs/1e9 0 1.2*max(abs(eq.Hrx).^2)])
end
    

