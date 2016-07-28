function [H, W] = calculate_fiber_preemphasis(N, lamb, Fiber, Hrx, sim)
    c = 299792458;
    tol = 1e-2;

%     Dtotal = 0;
%     for k = 1:length(Fibers)
%         Dtotal = Dtotal + Fibers(k).D(lamb)*Fibers(k).L;
%     end
%         
%     if Dtotal ~= 0  
%         theta = 1/2*((lamb.^2)/(2*pi*c))*(2*pi*sim.f).^2*Dtotal; % theta = -1/2*beta2*w.^2*L
%         Himdd = cos(theta);  
        Himdd = Fiber.Himdd(sim.f, lamb, 0, 'small signal');
        
        H = Himdd.*Hrx;
        h = fftshift(real(ifft(ifftshift(H))));
        s = cumtrapz(sim.t, abs(h).^2)/trapz(sim.t, abs(h).^2);
        htrim = h;
        htrim(s < tol & s > 1-tol) = [];
        h = htrim(1:sim.Mct/sim.ros.txDSP:end);
        
%         H = freqz(h, 1, ff, fs);
        fs = sim.ros.txDSP*sim.fs/sim.Mct;
        f = freq_time(length(h), fs);
        
%         H = fftshift(fft(h));
        H = freqz(h, 1, f, fs).*exp(1j*2*pi*f/fs*(length(h)-1)/2);
        H = H/sqrt(interp1(f, abs(H).^2, 0));
        H = 1./H;
        
        W = firls(N, f(f >= 0)/(fs/2), abs(H(f >=0 )));        
%         Hw = freqz(W, 1, f, fs);
        
        figure, hold on
        plot(f/1e9, 1./abs(H).^2)
%         plot(f/1e9, 1./abs(Hw).^2)  
        drawnow
end
    
    