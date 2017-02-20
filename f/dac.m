function [xa, xqzoh, xq] = dac(x, DAC, sim, verbose)
%% Digital-to-analog conversion: quantize, sample and hold with ZOH, and analog filtering
% Inputs:
% - x : input signal already at the DAC sampling rate
% - DAC : struct with DAC parameters. 
%    - ros : oversampling ratio of DAC with respect to symbol rate
%    - resolution : effective DAC resolution in bits. Quantization is only
%    applied if sim.quantiz = true and DAC.resolution not Inf
%    - filt : DAC filter after ZOH. Designed with design_filter.m. If filt
%    is empty then, filtering is not performed.
%    - offset (optional, default = 0) : additional time offset (in number 
%    of samples) to be removed i.e., xa(t) = xa(t - offset/fs). 
%       * Offset is only removed if filt is not empty.
%       * Offset need not be integer.
%    - rclip (optional, default = 0): clipping ratio. Clipping ratio is 
%    defined as the percentage of the signal amplitude that is clipped in 
%    both extremes (see code).
%    - excursion (optional, default = []): excursion limits of DAC. Specify
%    as a [xmin, xmax]. 
% - sim : simulation parameters
%     - f : frequency vector
%     - fs : sampling frequency to emulate continuous time
%     - Mct : oversampling ratio to emulate continuous time
% - verbose (optional, default = false) : whether to plot eye diagram
% Outputs:
% - xa : output signal at rate sim.fs
% - xqzoh : signal before analog filtering (after ZOH)

% Quantize
if isfield(sim, 'quantiz') && sim.quantiz && ~isinf(DAC.resolution)
    enob = DAC.resolution;
    if isfield(DAC, 'rclip')
        rclip = DAC.rclip;
    else
        rclip = 0;
    end
    
    if isfield(DAC, 'excursion') && not(isempty(DAC.excursion))
        xmax = DAC.excursion(2);
        xmin = DAC.excursion(1);
    else
        xmax = max(x);
        xmin = min(x);
    end
    xamp = xmax - xmin;    
    
    % Clipping
    % clipped: [xmin, xmin + xamp*rclip) and (xmax - xamp*rclip, xmax]
    % not clipped: [xmin + xamp*rclip, xmax - xamp*rclip]
    xmin = xmin + xamp*rclip; % discounts portion to be clipped
    xmax = xmax - xamp*rclip;
    xamp = xmax - xmin;
    
    dx = xamp/(2^(enob)-1);
    Pc = mean(x > xmax | x < xmin);
    if Pc ~= 0
        fprintf('DAC: clipping probability = %G\n', Pc) 
    end
    
    codebook = xmin:dx:xmax;
    partition = codebook(1:end-1) + dx/2;
    [~, xq, varQ] = quantiz(x, partition, codebook); 
    fprintf('DAC: quantization noise variance = %G | Signal-to-quantization noise ratio = %.2f dB\n', varQ, 10*log10(var(x)/varQ));
else
    xq = x;
    varQ = 0;
end

% Zero-order holder (ZOH)
Nhold = sim.Mct/DAC.ros; % number of samples to hold
assert(floor(Nhold) == ceil(Nhold), 'dac: oversampling ratio of DAC (DAC.ros) must be an integer multiple of oversampling ratio of continuous time (sim.Mct)');
xqzoh = upsample(xq, Nhold);
xqzoh = filter(ones(1, Nhold), 1, xqzoh); % group delay due to ZOH is removed inn the frequency domain

% Filtering
if not(isempty(DAC.filt))
    % Filter
    Hdac = ifftshift(DAC.filt.H(sim.f/sim.fs)); % rx filter frequency response
    
    % Time shift. Removes any group delay specified by offset. This is
    % equivalent to resampling the signal
    if isfield(DAC, 'offset')
        offset = DAC.offset + (Nhold-1)/2; % time offset in number of samples. Need not be integer
    else
        offset = (Nhold-1)/2; % if offtset was not defined only removed offset due to ZOH
    end
    
    Hshift = ifftshift(exp(1j*2*pi*sim.f/sim.fs*offset));

    % Filtering
    xa = real(ifft(fft(xqzoh).*Hdac.*Hshift)); % filter
else
    xa = xqzoh;
end
                                
% Plot
if exist('verbose', 'var') && verbose
    Ntraces = 100;
    Nstart = sim.Ndiscard*sim.Mct+1;
    Nend = min(Nstart + Ntraces*2*sim.Mct, length(xa));
    figure(301)
    subplot(221), box on
    eyediagram(xqzoh(Nstart:Nend), 2*sim.Mct)
    title('Eye diagram after ZOH')
    subplot(222), box on
    eyediagram(xa(Nstart:Nend), 2*sim.Mct)
%     eyediagram(xa(Nstart:Nend), 2*sim.Mct, sim.Mct, ceil(sim.Mct/2))
    title('DAC output eye diagram')
    subplot(223), box on
    plot(sim.f/1e9, abs(fftshift(fft(xa-mean(xa)))).^2)
    xlabel('Frequency (GHz')
    ylabel('|X(f)|^2')
    title('DAC output spectrum')
    a = axis;
    axis([-DAC.fs/1e9 DAC.fs/1e9 a(3:4)])
    drawnow
end