function pulse_shape = select_pulse_shape(type, sps, varargin)
%% Generate struct pulse_shape containing pulse shaping parameters
% Inputs:
% type: either 'rect' (rectangular), 'rrc' (root-raised cosine), or 'rc'
% (raised cosine)
% sps: samples per symbol
% rolloff: rolloff factor of pulse shape (0, 1]. Only required if type = 'rrc' or 'rc'
% span: number symbols over which pulse shape spans. Only required if type = 'rrc' or 'rc'

pulse_shape.type = type; % either 'rect', 'rrc', or 'rc'
pulse_shape.sps = sps; % number of samples per symbol
if not(strcmpi(type, 'rect'))
    assert(length(varargin) == 2, 'select_pulse_shape: invalid number of inputs. If type = rc or rrc, rolloff and span must be provided.')
end

switch lower(pulse_shape.type) 
    case 'rect' % Rectangular pulse shape
        pulse_shape.rolloff = 1;
        pulse_shape.h = ones(1, pulse_shape.sps); % pulse shapping filter coefficients
    case 'rrc' % Root raised cosine pulse shape
        pulse_shape.rolloff = varargin{1}; % 
        pulse_shape.span = varargin{2}; % number of symbols over which pulse shape spans
        pulse_shape.h = rcosdesign(pulse_shape.rolloff, pulse_shape.span, pulse_shape.sps, 'sqrt'); % pulse shapping filter coefficients
    case 'rc' % Raised cosine pulse shape
        pulse_shape.rolloff = varargin{1}; % 
        pulse_shape.span = varargin{2}; % number of symbols over which pulse shape spans
        pulse_shape.h = rcosdesign(pulse_shape.rolloff, pulse_shape.span, pulse_shape.sps, 'normal'); % pulse shapping filter coefficients
    otherwise
        error('select_pulse_shape: Invalid pulse shape type')
end