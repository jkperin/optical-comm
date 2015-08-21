%% Design equalizer and equalize yt
function [yd, eq] = equalize(type, yt, mpam, tx, fiber, rx, sim)
% Types supported:
% - None: No equalizer. Filter using antialiasing filter (rx.elefilt) and
% sample at symbol rate.
% - Analog: Invert channel response and apply matched-filter matched to the
% pulse shape and then do symbol-rate sampling
% - Fixed or Adaptive TD-SR-LE: time-domain symbol-rate linear equalizer.
% This assumes that there is a matched filter before symbol-rate sampling.
% - Fixed or Adaptive TD-FS-LE: time-domain fractionally-spaced linear equalizer.
% This assumes an anti-aliasing filter before sampling.

eq = rx.eq;

% Check if number of taps was defined and is even. If even, make it odd
if isfield(eq, 'Ntaps') && mod(eq.Ntaps, 2) == 0
    eq.Ntaps = eq.Ntaps + 1;
end

% If modulator wasn't defined in tx, don't do equalization
if ~isfield(tx, 'modulator')
    type = 'None';
end

% Inserts delay of half of sample if oversampling of continuous time is even
if mod(sim.Mct, 2) == 0 
    Delay = exp(-1j*pi*sim.f/sim.fs);
else 
    Delay = 1;
end

switch type
    case 'None'
        %% No equalization
        eq.H = ones(size(sim.f));
        eq.f = sim.f;
               
        % Reciver filter
        yt = real(ifft(fft(yt).*ifftshift(Delay.*rx.elefilt.H(sim.f/sim.fs))));
        
        yd = yt(floor(sim.Mct/2)+1:sim.Mct:end); % +1 because indexing starts at 1
        
        eq.Kne = 1;
        
    case 'Analog'
        %% Analog Equalization
        Gtx = design_filter('matched', mpam.pshape, 1/sim.Mct); % transmitter frequency response
        
        % Channel response
        Gch = tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay)...
            .*fiber.Hfiber(sim.f, tx);
        
        % Antialiasing filter
        Haa = design_filter('bessel', 5, 1/2*mpam.Rs/(sim.fs/2));
        Haa = Haa.H(sim.f/sim.fs);
               
        %        
        Heq = Haa;
        Heq(abs(sim.f)<=mpam.Rs/2) = Heq(abs(sim.f)<=mpam.Rs/2)./(Gch(abs(sim.f)<=mpam.Rs/2));

        % Apply equalization filter
        if isempty(yt) % function was called just to calculate equalizer
            yd = [];
        else
            % Channel inversion
            yteq = real(ifft(fft(yt).*ifftshift(Heq))); 
            
            % Matched filter
            yteq = real(ifft(fft(yteq).*ifftshift(Delay.*conj(Gtx.H(sim.f/sim.fs)))));

            % Sample
            yd = yteq(floor(sim.Mct/2)+1:sim.Mct:end); % +1 because indexing starts at 1
        end
        
        eq.H = Heq;
        eq.f = sim.f;
        
        eq.Kne = trapz(sim.f, abs(Heq).^2)/(mpam.Rs);
        
    case 'Fixed TD-FS-LE'
%             Grx = rx.elefilt;
            
            % Channel frequency response
%             Hch = tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay)...
%                 .*fiber.Hfiber(sim.f, tx).*Grx(sim.f/sim.fs);
%             
% %             yt = real(ifft(fft(yt).*ifftshift(Hmatched)));
% % 
% %             yk = yt(floor(sim.Mct/2)+1:sim.Mct:end);
%             
%             %% MMSE Time-domain symbol-rate equalizer
%             gtx = mpam.pshape(sim.t*sim.fs + (sim.Mct-1)/2);
%             hch = ifft(ifftshift(Hch)); % received pulse shape
%             p = conv(gtx, hch, 'same');
%             p = p/max(abs(p));
%                         
%             
%             plot(sim.t, p)
%             1;
                        
%             xd = x(mod(abs(n), sim.Mct) == 0); % symbol-rate sample received pulse
%             
%             
%             xd = xd/abs(sum(xd)); % normalize to unit gain at DC
%             
%             X = toeplitz([xd.'; zeros(eq.Ntaps-1, 1)], [xd(1) zeros(1, eq.Ntaps-1)]);
%             
%             % Get correct number of taps from Toeplitz matrix
%             X = X(ceil(size(X, 1)/2)+(-floor(eq.Ntaps/2):floor(eq.Ntaps/2)), :);
%             
%             e = zeros(eq.Ntaps, 1); 
%             e((eq.Ntaps-1)/2) = 1;
%             
%             % NSR = noise signal ratio
%             if isfield(eq, 'NSR')
%                 W = X*(((X' + eq.NSR*eye(eq.Ntaps))*X)\e);
%             else
%                 W = X*((X'*X)\e);
%             end                
%             % Note: if NSR = 0, or if NSR is not specified, then
%             % zero-forcing equalizer is designed            
% 
%             % Filter using
%             if isempty(yt) % only design
%                 yd = []; 
%             else
%                 % Matched filter
%                 yt = real(ifft(fft(yt).*ifftshift(Hmatched)));
%                 % Sample at symbol rate
%                 yk = yt(floor(sim.Mct/2)+1:sim.Mct:end);
%                 % Filter
%                 yd = filter(W, 1, yk);
%                 yd = circshift(yd, [-(eq.Ntaps-1)/2+1 0]); % remove delay of transversal filter
%             end  
% 
%             % Aux
%             eq.num = W;
%             eq.den = 1;
%             [eq.H, w] = freqz(eq.num, eq.den);
%             eq.f = w/(2*pi);
%             eq.Kne = 2*trapz(eq.f, abs(eq.H));
        
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
        eq.Kne = trapz(eq.f, abs(eq.H).^2);
        
    case 'Adaptive TD-SR-LE'
        Ntaps = eq.Ntaps;
        Ntrain = eq.Ntrain;
        mu = eq.mu;
        b = mpam.mod(rx.eq.TrainSeq, 1);
        
        %% Time-domain symbol-rate equalizer
        % matchedfilt = design_filter('matched', @(t) conv(mpam.pshape(t), 1/sim.fs*tx.modulator.h(t/sim.fs), 'full') , 1/sim.Mct); 
        Gtx = design_filter('matched', mpam.pshape, 1/sim.Mct); 

        % Received Pulse Spectrum
        Hmatched = Delay.*conj(Gtx.H(sim.f/sim.fs).*tx.modulator.H(sim.f).*...
            exp(1j*2*pi*sim.f*tx.modulator.grpdelay).*fiber.Hfiber(sim.f, tx));

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
            Gtx = design_filter('matched', mpam.pshape, 1/sim.Mct); 

            Hmatched = Delay.*conj(Gtx.H(sim.f/sim.fs).*tx.modulator.H(sim.f).*...
                exp(1j*2*pi*sim.f*tx.modulator.grpdelay).*fiber.Hfiber(sim.f, tx));
            
%             yt = real(ifft(fft(yt).*ifftshift(Hmatched)));
% 
%             yk = yt(floor(sim.Mct/2)+1:sim.Mct:end);
            
            %% MMSE Time-domain symbol-rate equalizer
            n = -floor(eq.Ntaps/2)*sim.Mct:sim.Mct*floor(eq.Ntaps/2);
            v = mpam.pshape(n);
            h = tx.modulator.h(n/sim.fs + tx.modulator.grpdelay);
            p = conv(h, v, 'same');
            x = conv(p, fliplr(p), 'same');
            xd = x(mod(abs(n), sim.Mct) == 0); % symbol-rate sample received pulse
            xd = xd/abs(sum(xd)); % normalize to unit gain at DC
            
            X = toeplitz([xd.'; zeros(eq.Ntaps-1, 1)], [xd(1) zeros(1, eq.Ntaps-1)]);
            
            % Get correct number of points from Toeplitz matrix
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
            
%             % get filter coefficients
%             W = conj(W(end:-1:1));

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
            
            %% Design of zero-forcing TD-SR_L             
%             % Fold spectrum
%             df = sim.fs/length(sim.f);
%             X = abs(Hmatched).^2;
%             X = conj(flipud(X(1:floor(length(X)/2)+1)));
%             ff = -flipud(sim.f(1:floor(length(sim.f)/2)+1));
%             Nfold = floor(mpam.Rs/(2*df))+1;
%             Xfolded = X(1:Nfold);
%             ffolded = ff(1:Nfold);
%             fold = true;
%             for k = Nfold+1:Nfold:length(ff)-Nfold
%                 if fold
% %                     [ff(k+Nfold-1) ff(k)]/1e9
% 
%                     Xfolded = Xfolded + conj(flipud(X(k:k+Nfold-1)));
% 
%                     fold = false;
%                 else
% %                     [ff(k) ff(k+Nfold-1)]/1e9
% 
%                     Xfolded = Xfolded + X(k:k+Nfold-1);
% 
%                     fold = true;
%                 end
%             end
% 
%             Xfolded = [flipud(conj(Xfolded(2:end))); Xfolded(1:end-1)];
%             ffolded = [-flipud((ffolded(2:end))); ffolded(1:end-1)];
%            
%             dff = 1/length(yk);
%             ff = -0.5:dff:0.5-dff;
%             ff = mpam.Rs*ff.';
%             
%             eq.f = ff;
%             eq.H = interp1(ffolded, 1./Xfolded, ff, 'spline');
%             
%             % filter using Xfolded
%             yd = real(ifft(fft(yk).*ifftshift(eq.H)));
%             
%             eq.Kne = trapz(eq.f, abs(eq.H))/(mpam.Rs);
        
    otherwise
        error('Equalization option not implemented yet!')
end

if sim.verbose
    plot(sim.f/1e9, abs(Heq).^2)
    xlabel('Frequency (GHz)')
    ylabel('|H_{eq}(f)|^2')
    axis([0 2*mpam.Rs/1e9 0 1.2*max(abs(Heq).^2)])
end
    

