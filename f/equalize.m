function [yd, Heq] = equalize(type, yt, mpam, tx, fiber, rx, eq, sim)

switch type
    case 'None'
        %% No equalization
        Heq = ones(size(sim.f));
               
        % Antialiasing filter
        yt = real(ifft(fft(yt).*ifftshift(rx.elefilt.H(sim.f/sim.fs))));
        
        yd = yt(floor(sim.Mct/2):sim.Mct:end);
        
    case 'Analog'
        %% Analog Equalization
        Gtx = design_filter('matched', mpam.pshape, 1/sim.Mct); % transmitter frequency response
        
        % Received Pulse Spectrum
        Pf = Gtx.H(sim.f/sim.fs).*tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay)...
            .*fiber.Hfiber(sim.f, tx);
        
        % Antialiasing filter
        Haa = design_filter('bessel', 5, 1/2*mpam.Rs/(sim.fs/2));
%         Haa = rx.elefilt.H(sim.f/sim.fs);
        Haa = Haa.H(sim.f/sim.fs);
        
        %
%         Pf = Pf(abs(sim.f) <= mpam.Rs/2); % take Pf response that falls inside the max frequency
%         Heq = Haa;
%         Heq(abs(sim.f) <= mpam.Rs/2) = Heq(abs(sim.f) <= mpam.Rs/2)./Pf;
%         Heq(abs(sim.f) >  mpam.Rs/2) = Pf(end)*Haa(abs(sim.f) > mpam.Rs/2);
        
        %        
        Heq = Haa;
        Heq(abs(sim.f)<=mpam.Rs/2) = Heq(abs(sim.f)<=mpam.Rs/2)./(Pf(abs(sim.f)<=mpam.Rs/2));

        % Apply equalization filter
        if isempty(yt) % function was called just to calculate equalizer
            yd = [];
        else
            yteq = real(ifft(fft(yt).*ifftshift(Heq)));

            % Sample
            yd = yteq(floor(sim.Mct/2):sim.Mct:end);
        end
        
    case 'Fixed TD-FSE'
    case 'Adaptive TD-FS-LE'
        %% Adaptive Time-domain fractionally-spaced equalizer
        ros = eq.ros;
        Ntaps = eq.Ntaps;
        Ntrain = eq.Ntrain;
        mu = eq.mu;
        b = mpam.mod(eq.TrainSeq, 1);
        
        % Antialiasing filter
        ytaa = real(ifft(fft(yt).*ifftshift(rx.elefilt.H(sim.f/sim.fs))));
        
        if mod(sim.Mct/ros, 2) == 0
            yk = ytaa(1:sim.Mct/ros:end);
            tk = sim.t(1:sim.Mct/ros:end);
        else
%             yk = interp1(1:length(ytaa), ytaa, 1:sim.Mct/ros:length(ytaa));
%             tk = interp1(1:length(ytaa), sim.t, 1:sim.Mct/ros:length(ytaa));
            
            yk = resample(ytaa, sim.ros, sim.Mct);
            tk = resample(sim.t, sim.ros, sim.Mct);            
        end
        
        W = [zeros(Ntaps-1, 1); 1]; % Filter taps
        y = zeros(size(yk));
        n = 1;
        e = zeros(1, length(yk)/ros);
        for k = max(Ntaps, sim.Ndiscard*ros):length(yk)
            z = yk(k-Ntaps+1:k);

            y(k) = sum(W.*z);

            if mod(k, ros) == 0
                if n < Ntrain % Training
                    e(k/ros) = y(k) - b(k/ros);
                else
                    e(k/ros) = y(k) - mpam.mod(mpam.demod(y(k)), 1);
                end
                n = n + 1;

                W = W - 2*mu*e(k/ros)*z;
            end
        end

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
        
        Heq = W; % return taps instead of frequency response
    otherwise
        error('Equalization option not implemented yet!')
end

if sim.verbose
    plot(sim.f/1e9, abs(Heq).^2)
    xlabel('Frequency (GHz)')
    ylabel('|H_{eq}(f)|^2')
    axis([0 2*mpam.Rs/1e9 0 1.2*max(abs(Heq).^2)])
end
    

