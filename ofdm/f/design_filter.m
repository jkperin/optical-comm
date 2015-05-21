% input
% type = filter type {'ideal', 'butter', 'cheby1', 'ellipt', 'gaussian',
% 'bessel', 'fir'}
% order = filter order
% cutoff = cutoff frequency normalized with respect to the sampling rate of
% the ofdm signal fs (ofdm.fs)
% Mct = oversampling ratio to simulate continuous time

% output
% filt = impulse response of the filter
% bfilt = numerator coefficients of filter transfer functions
% afilt = denominator coefficients of filter transfer functions

function [bfilt, afilt] = design_filter(type, order, cutoff, Mct)
switch type
    case 'ideal'
        % ideal interpolation using interp
        afilt = 1;
        bfilt = intfilt(Mct, 9, cutoff);  
        
    case 'butter'
        % Butterworth filter
        [bfilt, afilt] = butter(order, cutoff/Mct);
        
    case 'cheby1'
        % Chebyshev 1 filter
        ripple = 0.1;
        
        [bfilt, afilt] = cheby1(order, ripple, cutoff/Mct);
        
    case 'ellipt'
        % Elliptic 1 filter
        ripple = 0.1;
        stopband_att = 50;
        
        [bfilt, afilt] = ellip(order, ripple, stopband_att, cutoff/Mct);
        
    case 'gaussian'
        nlength = 256;   % number of samples used for fitting. Verified empirically

        % Gaussian filter
        df = 1/nlength;
        f = -0.5:df:0.5-df;
        f0 = (cutoff/Mct)/(2*log(2)^(1/(2*order)));
        H = exp(-0.5*(f/f0).^(2*order));

        h = real(fftshift(ifft(ifftshift(H))));

        afilt = 1;
        bfilt = h;
        
%         freqz(bfilt, afilt, f, 1)
%         1;
%        
    case 'bessel'       
        % Generates CT Bessel filter prototype and convert it to DT using
        % the bilinear transformation with frequency prewarping.
        % Note 1: this is equivalent to what matlab does for the other
        % filters (i.e., DT butter is equal to taking bilinear transformation with frequency prewarping 
        % of the butter filter designed in CT).
        % Note 2: The wo parameter in Bessel filter design is the frequency 
        % up to which the filter's group delay is approximately constant.
        % Factor of wc2wo converts cutoff frequency into wo. 
        wc2wo = 1.621597623980423; 
        % Note: wc2w0 was calculated using the following procedure:
        % function h = calc_besself_cutoff(order, f3dB, wc2wo)
        % 
        %     [b,a] = besself(order, 2*pi*f3dB*abs(wc2wo));
        % 
        %     h = polyval(b, 1j*2*pi*f3dB)/polyval(a, 1j*2*pi*f3dB);
        % 
        %     h = abs(h).^2;
        % 
        % Calculate wc2w0 that makes h == 0.5
        % fzero(@(x) calc_besself_cutoff(5, 1, x) - 0.5, 1.65)
        % The value of x doens't depend on f3dB

        wow = 2*tan(pi*wc2wo*cutoff/(2*Mct)); % prewarped cutoff frequency
        [bfilt, afilt] = besself(order, wow); % design prototype in CT with frequency prewarped
        [bfilt, afilt] = bilinear(bfilt, afilt, 1);
        
%         sys = tf(bfilt, afilt);
%         sysz = c2d(sys, 1, 'tustin'); % convert into DT by bilinear transformation
%         [bfilt, afilt] = tfdata(sysz);
%         bfilt = cell2mat(bfilt);
%         afilt = cell2mat(afilt);  
        
    case 'fir'
        % FIR
        imp_length = 100;
                
        bfilt = fir1(imp_length-Mct, cutoff/Mct);
        afilt = 1;
        freqz(bfilt, afilt);
       
    otherwise
        assert(false, 'Unknown option!')
end

