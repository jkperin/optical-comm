function dlamb = df2dlamb(df, range)
%% Convert variation in frequency (df) to variation in wavelength (dlamb) near range (range)
% Inputs:
% - df: variation in frequency in Hz
% - range (optional, default = 1550): either 1310 or 1550
% Output:
% - dlamb: variation in wavelength near range
if not(exist('range', 'var'))
    range = 1550;
end

c = 299792458;      % speed of light
range = range*1e-9; % convert to nm
f1 = c/range;
f2 = f1 + df;
lamb2 = c/f2;
dlamb = 1550e-9-lamb2;
