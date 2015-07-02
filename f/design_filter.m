%% Design filters
% Continuous-time filters are designed in continuous time and then
% converted to discrete time using bilinear transformation with frequency
% prewarping. If oversampling ratio of 'continuous time' (Mct) is high enough
% approximation will be accurate.

% Gaussian and FBG filters are converted to discrete time by impulse
% invariance

% input
% type = filter type {'butter', 'cheby1', 'ellipt', 'gaussian',
% 'bessel', 'fir', 'lorentzian', 'fbg'}
% order = filter order
% fcnorm = normalized cutoff frequency. fcnorm = f3dB/(fs/2). 0 < fcnorm < 1
% Mct = oversampling ratio to simulate continuous time

% output
% struct filt with
% a : denominator coefficients
% b : numerator coefficients
% grpdelay : group delay
% H : anonymous function of frequency response with group delay eliminated
% noisebw : anonymous function of equivalent two-sided noise bandwidth over
% larger number of samples 2^15 (large number) given a sampling frequency fs 

function filt = design_filter(type, order, fcnorm)
switch type       
    case 'butter'
        % Butterworth filter
        [num, den] = butter(order, fcnorm);
        
    case 'cheby1'
        % Chebyshev 1 filter
        ripple = 0.1;
        
        [num, den] = cheby1(order, ripple, fcnorm);
        
    case 'ellipt'
        % Elliptic 1 filter
        ripple = 0.1;
        stopband_att = 50;
        
        [num, den] = ellip(order, ripple, stopband_att, fcnorm);
        
    case 'matched'
        % order = pulse shape function
        % fcnorm = 1 / oversampling factor
        nmax = 1024; % more than enough 
        pshape = order;
        cumenergy = cumsum(abs(pshape(0:nmax-1)).^2)/sum(abs(pshape(0:nmax-1)).^2);
        nlength = find(cumenergy >= 0.9999999, 1, 'first');
        
        num = conj(fliplr(pshape(0:nlength-1)));
        num = num/sum(num); % ensures unit gain at DC
        den = 1;
    
    case 'fir'
        % FIR
        imp_length = 100;
                
        num = fir1(imp_length, fcnorm);
        den = 1;
        freqz(num, den); 
        
    case 'gaussian'
        nlength = 256;   % number of samples used for fitting. 

        % Gaussian filter
        df = 1/nlength;
        f = -0.5:df:0.5-df;
        f0 = (fcnorm)/(2*log(2)^(1/(2*order)));
        H = exp(-0.5*(f/f0).^(2*order));

        h = real(fftshift(ifft(ifftshift(H))));

        den = 1;
        num = h;
        
    case 'bessel'       
        % Generates CT Bessel filter prototype and convert it to DT using
        % the bilinear transformation with frequency prewarping.
        % Note 1: this is equivalent to what matlab does for the other
        % filters (i.e., DT butter is equal to taking bilinear transformation with frequency prewarping 
        % of the butter filter designed in CT).
        % Note 2: The wo parameter in Bessel filter design is the frequency 
        % up to which the filter's group delay is approximately constant.
        % Factor of wc2wo converts fcnorm frequency into wo. 
        wc2wo = 1.621597623980423; 
        % Note: wc2w0 was calculated using the following procedure:
        % function h = calc_besself_fcnorm(order, f3dB, wc2wo)
        % 
        %     [b,a] = besself(order, 2*pi*f3dB*abs(wc2wo));
        % 
        %     h = polyval(b, 1j*2*pi*f3dB)/polyval(a, 1j*2*pi*f3dB);
        % 
        %     h = abs(h).^2;
        % 
        % Calculate wc2w0 that makes h == 0.5
        % fzero(@(x) calc_besself_fcnorm(5, 1, x) - 0.5, 1.65)
        % The value of x doens't depend on f3dB

        wow = wc2wo*(2*pi*fcnorm);
        [nums, dens] = besself(order, wow); % design prototype in CT with frequency prewarped   
%         [h, w]= freqs(num, den, 2^10);
%         plot(w/(2*pi), abs(h).^2)
        
        [num, den] = bilinear(nums, dens, 2, fcnorm);
%         [h, w] = freqz(num, den);
%         hold on, 
%         plot(w/(2*pi), abs(h).^2);
%         1;

                        
    case 'lorentzian'
        den = [2/(2*pi*2*fcnorm) 1];
        num = 1;
        [num, den] = bilinear(num, den, 2, fcnorm);

    case 'fbg' % fiber Bragg grating
        % Values based on http://ee.stanford.edu/~jmk/pubs/dpsk.cd.pmd.pdf
        kappa = 6; % cm^-1
        Lg = 2/kappa; %
        betaf = @(f, vg) sqrt((2*pi*f/vg).^2 - kappa^2);
        Hf = @(f, vg) 1/tanh(kappa*Lg)*1j*kappa.*sin(betaf(f, vg)*Lg)./(betaf(f, vg).*cos(betaf(f, vg)*Lg) - 1j*2*pi*(f/vg).*sin(betaf(f, vg)*Lg));
        
        [vg, ~, exitflag] = fzero(@(vg) abs(Hf(fcnorm/2, abs(vg))).^2 - 0.5, fcnorm);
        vg = abs(vg);
        
        if exitflag ~= 1
            warning('fiber Bragg grating filter did not converge to desired bandwidth')
        end
        
        N = 2^8;
        df = 1/N;
        f = -0.5:df:0.5-df; 
               
        den = 1;
        num = fftshift(ifft(ifftshift(Hf(f, vg))));  
           
    otherwise
        assert(false, 'Unknown option!')
end

filt.num = num;
filt.den = den;
filt.grpdelay = grpdelay(num, den, 1);
filt.fcnorm = fcnorm;
filt.H = @(f) freqz(num, den, 2*pi*f).*exp(1j*2*pi*f*filt.grpdelay);
filt.noisebw = @(fs) noisebw(num, den, 2^15, fs); % equivalent two-sided noise bandwidth over larger number of samples (2^15) given a sampling frequency fs 


% figure, 
% subplot(211)
% plot(f, abs(filt.H(f).^2))
% subplot(212)
% plot(f, unwrap(angle(filt.H(f))))