%% Design filters. See validate_design_filter.m
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
% num : denominator coefficients
% den : numerator coefficients
% grpdelay : group delay
% H : anonymous function of frequency response with group delay eliminated
% h : vector containing impulse response with unit gain at DC
% noisebw : anonymous function of equivalent two-sided noise bandwidth over
% larger number of samples 2^15 (large number) given a sampling frequency fs 

function filt = design_filter(type, order, fcnorm, verbose)

if not(exist('verbose', 'var'))
    verbose = false;
end

% Variables used when calculating impulse response
maxMemoryLength = 2^10; % maximum memory length
threshold = 1-1e-6;      % energy threshold above which impulse respone is neglected

type = lower(type); % converts to lowercase
switch type       
    case 'butter'
        % Butterworth filter
        [num, den] = butter(order, fcnorm);
            
    case 'cheby1'
        % Chebyshev 1 filter
        ripple = 0.5;
        
        [num, den] = cheby1(order, ripple, fcnorm);
        
    case 'ellipt'
        % Elliptic 1 filter
        ripple = 0.1;
        stopband_att = 50;
        
        [num, den] = ellip(order, ripple, stopband_att, fcnorm);
        
    case 'two-pole'
        % 2nd-order filter with unit damping
        % The continuous-time transfer function is converted to
        % discrete-time using the bilinear transformation with frequency
        % prewarping
        bw = order; % the argument order is used as bandwidth
        fs = fcnorm; % the argument fcnorm is used as sampling frequency
        order = 2;
        fcnorm = bw/(fs/2);
        
        wc = 2*pi*bw/sqrt(sqrt(2)-1); % converts BW into fc
        bs = 1;
        as = [(1/wc).^2, 2/wc, 1];
        [num, den] = bilinear(bs, as, fs, bw);
 
        if verbose
            verbose = false;
            f = linspace(0, fs/2);
            Hs = freqs(bs, as, 2*pi*f);
            Hz = freqz(num, den, f, fs);
            plot_transform(Hs, Hz, f/fs, bw/fs);
        end

    case 'matched'
        % order = pulse shape function
        % fcnorm = 1 / oversampling factor
        nmax = 256; % more than enough 
        pshape = order;
        Mct = 1/fcnorm;
        
        cumenergy = cumsum(abs(pshape(0:nmax-1)).^2)/sum(abs(pshape(0:nmax-1)).^2);
        nlength = find(cumenergy >= 0.9999, 1, 'first');
        
        num = conj(pshape(nlength-1:-1:0));
        num = num/sum(num); % ensures unit gain at DC
        den = 1;
        
        % Note: if filtering with filter(num, den, x) the delay inserted by
        % the transversal filter implementation must be removed, otherwise
        % the matched filter will maximize the SNR not at the sampling
        % instant (i.e., middle of the pulse)
        
        % Noise bandwidth
        nbw = @(fs) sum(abs(num).^2)*fs; % requires unit gain at DC
        
    case 'fir1'
        % FIR
        imp_length = order;
                
        num = fir1(imp_length, fcnorm);
        den = 1;
        
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
        
        % Noise bandwidth
        nbw = @(fs) sum(abs(num).^2)*fs; % requires unit gain at DC
        
        if verbose
            verbose = false;
            f = linspace(0, 1/2);
            Hs = exp(-0.5*(f/f0).^(2*order));
            Hz = freqz(num, den, f, 1);
            plot_approx(Hs, Hz, f, fcnorm/2)
        end
                
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

        wow = wc2wo*(2*pi*fcnorm/2);
        [nums, dens] = besself(order, wow); % design prototype in CT with frequency prewarped           
%         [num, den] = bilinear(nums, dens, 1, fcnorm/2);
        [num, den] = impinvar(nums, dens, 1);

        if verbose
            verbose = false;
            f = linspace(0, 1/2);
            Hs = freqs(nums, dens, 2*pi*f);
            Hz = freqz(num, den, f, 1);
            plot_transform(Hs, Hz, f, fcnorm/2)
        end

    case 'lorentzian'
        dens = [2/(2*pi*2*fcnorm/2) 1];
        nums = 1;
        [num, den] = bilinear(nums, dens, 1, fcnorm/2);
        
        if verbose
            verbose = false;
            f = linspace(0, 1/2);
            Hs = freqs(nums, dens, 2*pi*f);
            Hz = freqz(num, den, f, 1);
            plot_transform(Hs, Hz, f, fcnorm/2)
        end
        
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
        f = -1/2:df:1/2-df;
               
        den = 1;
        num = fftshift(ifft(ifftshift(Hf(f, vg)))); 
        
        % Noise bandwidth
        nbw = @(fs) sum(abs(num).^2)*fs; % requires unit gain at DC
        
        if verbose
            verbose = false;
            f = linspace(0, 1/2);
            Hs = Hf(f, vg);
            Hz = freqz(num, den, f, 1);
            plot_approx(Hs, Hz, f, fcnorm/2)
        end
           
    otherwise
        error('Unknown filter type = %s!', type)
end
filt.type = type;
filt.order = order;
filt.num = num;
filt.den = den;
filt.grpdelay = grpdelay(num, den, 1);
filt.fcnorm = fcnorm;
filt.H = @(f) freqz(num, den, 2*pi*f).*exp(1j*2*pi*f*filt.grpdelay);
if exist('nbw', 'var')
    filt.noisebw = nbw;
else
    filt.noisebw = @(fs) noisebw(num, den, 2^15, fs); % equivalent two-sided noise bandwidth over larger number of samples (2^15) given a sampling frequency fs 
end

if den == 1 % FIR
    filt.h = num;
else
    x = zeros(1, maxMemoryLength+1);
    x(1) = 1;
    y = filter(filt.num, filt.den, x);
    E = cumsum(abs(y).^2)/sum(abs(y).^2);
    y(E > threshold) = [];
    y = y/abs(sum(y)); % normalize to have unit gain at DC
    filt.h = y;
end

if verbose
    plot_filter(filt)
end
end

function plot_filter(filt)
    f = linspace(0, 0.5);
    Hz = filt.H(f);
    figure(111)
    subplot(211), hold on, box on
    plot(f, 20*log10(abs(Hz)), 'DisplayName', filt.type)
    aa = axis;
    h = plot([0 filt.fcnorm/2], [-3 -3], ':k');
    hasbehavior(h, 'legend', false);   % line will not be in legend
    h = plot(filt.fcnorm/2*[1 1], [aa(3) -3], ':k');
    hasbehavior(h, 'legend', false);   % line will not be in legend
    xlabel('Normalized frequency f/f_s')
    ylabel('Magnitude (dB)')
    legend('-DynamicLegend')
    axis([0 0.5 -50 10])
        
    subplot(212), hold on, box on
    plot(f, rad2deg(unwrap(angle(Hz))), 'DisplayName', filt.type)
    xlabel('Normalized frequency f/f_s')
    ylabel('Phase (deg)')
    legend('-DynamicLegend')
    axis([0 0.5 -360 360])
end

function plot_approx(Hs, Hz, f, fc)
    figure(112)
    subplot(211), hold on, box on
    plot(f, 20*log10(abs(Hs)), f, 20*log10(abs(Hz)), '--')
    aa = axis;
    plot([0 fc], [-3 -3], ':k')
    plot(fc*[1 1], [aa(3) -3], ':k')
    xlabel('Normalized frequency f/f_s')
    ylabel('Magnitude (dB)')
    legend('Analog filter', 'Approximation')
    axis([0 0.5 -50 10])

    subplot(212), hold on, box on
    plot(f, rad2deg(unwrap(angle(Hs))), f, rad2deg(unwrap(angle(Hz))), '--')
    xlabel('Normalized frequency f/f_s')
    ylabel('Phase (deg)')
    legend('Analog filter', 'FIR approximation')
    axis([0 0.5 -360 360])
end

function plot_transform(Hs, Hz, f, fc)
    figure(113)
    subplot(211), hold on, box on
    plot(f, 20*log10(abs(Hs)), f, 20*log10(abs(Hz)))
    aa = axis;
    plot([0 fc], [-3 -3], ':k')
    plot(fc*[1 1], [aa(3) -3], ':k')
    xlabel('Normalized frequency f/f_s')
    ylabel('Magnitude (dB)')
    legend('Analog filter', 'Bilinear transform')
    axis([0 0.5 -50 10])
    
    subplot(212), hold on, box on
    plot(f, rad2deg(unwrap(angle(Hs))), f, rad2deg(unwrap(angle(Hz))))
    xlabel('Normalized frequency f/f_s')
    ylabel('Phase (deg)')
    legend('Analog filter', 'Bilinear transform')
    axis([0 0.5 -360 360])
end
