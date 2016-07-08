function [ys, varQ, xa] = adc(x, ADC, sim, Hrx)
%% Analog-to-digital conversion: filter, downsample, and quantize signal x
% Inputs:
% - x : input signal
% - ADC : struct with ADC parameters. 
%    - ros : oversampling ratio of ADC with respect to symbol rate
%    - ENOB : effective number of bits. Quantization is only
%    applied if sim.quantiz = true and ADC.ENOB not Inf
%    - filt : ADC antialiasing filter. Designed with design_filter.m. If filt
%    is empty then, filtering is not performed.
%    - offset (optional, default = 0) : additional time offset (in number 
%    of samples) to be removed i.e., x(t) = x(t - offset/fs). 
%       * Offset is only removed if filt is not empty.
%       * Offset need not be integer.
%    - rclip (optional, default = 0): clipping ratio. Clipping ratio is 
%    defined as the percentage of the signal amplitude that is clipped in 
%    both extremes (see code).
% - sim : simulation parameters
%     - f : frequency vector
%     - fs : sampling frequency to emulate continuous time
%     - Mct : oversampling ratio to emulate continuous time
% - Hrx (optional) : if Hrx is provided ADC.filt is ignored and Hrx is
% applied isntead. Typically, Hrx is used when we want to do matched
% filtering rather than antialiasing filter.
% Outputs:
% - y : output signal
% - varQ : quantization noise variance
% - xa : signal after receiver filtering (before sampling)

% Filter
% Filter by lowpass filter (antialiasing filter)
if not(exist('Hrx', 'var')) % check if ADC filter is defined
    Hrx = ifftshift(ADC.filt.H(sim.f/sim.fs)); % rx filter frequency response
else % if Hrx exists use Hrx instead. This is typically a matched filter
    Hrx = ifftshift(Hrx); % time shift so that first sample correspond to the center of symbol
end

if isfield(ADC, 'offset') % time shift in number of samples. Need not be integer
    Hshift = ifftshift(exp(1j*2*pi*sim.f/sim.fs*ADC.offset)); % time shift to make first sample center of first pulse
else
    Hshift = 1;
end

% Filter
xa = real(ifft(fft(x).*Hrx.*Hshift));

% Downsample
if isInteger(sim.Mct/ADC.ros) % resampling is not required
    xs = xa(1:sim.Mct/ADC.ros:end);
else
    warning('adc: sim.Mct/ADC.ros is not interger, so resampling is performed')
    [N, D] = rat(ADC.ros);
    xs = resample(xa, D, N);
end

% Quantize
if isfield(sim, 'quantiz') && sim.quantiz && ~isinf(ADC.ENOB)
    enob = ADC.ENOB;
    if isfield(ADC, 'rclip')
        rclip = ADC.rclip;
    else
        rclip = 0;
    end
    
    xmax = max(xs);
    xmin = min(xs);
    xamp = xmax - xmin;    
    
    % Clipping
    % clipped: [xmin, xmin + xamp*rclip) and (xmax - xamp*rclip, xmax]
    % not clipped: [xmin + xamp*rclip, xmax - xamp*rclip]
    xmin = xmin + xamp*rclip; % discounts portion to be clipped
    xmax = xmax - xamp*rclip;
    xamp = xmax - xmin;
    
    dx = xamp/(2^(enob)-1);
    
    codebook = xmin:dx:xmax;
    partition = codebook(1:end-1) + dx/2;
    [~, ys, varQ] = quantiz(xs, partition, codebook); 
else
    ys = xs;
    varQ = 0;
end