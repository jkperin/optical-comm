function Nops = Nops_freq_domain_equalizer(Nfft, M, fft_type)
%% Number of operations per useful sample in frequency domain equalizer
% A frequency-domain implementation of a FIR filter of order M requires 2
% FFT operations and one multiplication and produces Nfft-M useful outputs.
% This applies to either overlap-save and overlap-add methods
% Inputs;
% - Nfft: FFT size
% - M: FIR filter order
% - fft_type(optional, default='split radix'): FFT implementation

if not(exist('fft_type', 'var'))
    fft_type = 'split radix';
end

switch lower(fft_type)
    case 'split radix'
        Offt = @(Nfft) 4*Nfft.*log2(Nfft); % number of real operations assuming split-radix FFT
    otherwise
        error('Nops_freq_domain_equalizer: invalid fft_type')
end

% Number of operations divided by number of useful samples
Nops = (2*Offt(Nfft) + Nfft)./(Nfft - M);