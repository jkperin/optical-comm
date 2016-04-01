function [ys, varQ] = adc(x, ADC, sim)
%% Digital to analog conversion: filters, downsample, and quantize signal x
% Inputs:
% - x : input signal
% - ADC : ADC parameters. Must contain fields {filt, ros}. 
% If sim.quantiz is true, then quantization is performed and ADC must contain
% the additional fields {ENOB, rclip}, ENOB and clipping ratio. 
% Clipping ratio is defined as the percentage of the signal amplitude that 
% is clipped in both extremes (see code)
% - sim : simulation parameters {sim.f, sim.fs, sim.Mct}
% Outputs:
% - y : output signal
% - varQ : quantization noise variance

% Filter
% Filter by lowpass filter (antialiasing filter)
Hrx = ifftshift(ADC.filt.H(sim.f/sim.fs).*... % rx filter frequency response
    exp(1j*2*pi*sim.f/sim.fs*(sim.Mct-1)/2)); % time shift so that first sample correspond to the center of symbol

xa = real(ifft(fft(x).*Hrx)); % filter                                   

% Downsample
xs = xa(1:sim.Mct/ADC.ros:end);

% Quantize
if isfield(sim, 'quantiz') && sim.quantiz && ~isinf(ADC.ENOB)
    enob = ADC.ENOB;
    rclip = ADC.rclip;
    
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


    

