%% Digital PLL analysis
clear, clc, close all

sim.Rs = 56e9;
Ts = 1/sim.Rs;

CT2DT = 'bilinear';                                                 % method for converting continuous-time loop filter to discrete time: {'bilinear', 'impinvar'}
Delay = 0;                                                          % (not implemented yet) Delay in DPLL loop measured in number of symbols 
csi = 1/sqrt(2);                                                    % damping coefficient of second-order loop filter
wn = 2*pi*0.8e9;        

if strcmpi(CT2DT, 'bilinear')
    wn = 2/Ts*atan(wn*Ts/2); % Compensates for frequency warping in bilinear transformation.
end

nums = [0 2*csi*wn wn^2];
dens = [1 0 0]; % descending powers of s
% Open-loop digital filter coefficients using bilinear transformation
if strcmpi(CT2DT, 'bilinear')
    [numz, denz] = bilinear(nums, dens, sim.Rs); % ascending powers of z^–1
elseif strcmpi(CT2DT, 'impinvar')
    [numz, denz] = impinvar(nums, dens, sim.Rs); % ascending powers of z^–1
else
    error('dpll/invalid method for continuous time to discrete time conversion')
end

% Continuous time
% Open loop frequency response
Gs = tf(nums, dens); 
Gz = tf(numz, denz, 1/sim.Rs);

%
figure, hold on
step(1/(Gs+1))
step(1/(Gz+1))
title('Phase error response to step in AWGN noise')
legend('ct', 'dt')

figure, hold on
step(Gs/(Gs+1))
step(Gz/(Gz+1))
title('Phase error response to step in phase noise')
legend('ct', 'dt')


figure, hold on
impulse(1/(Gs+1))
impulse(1/(Gz+1))
title('Phase error response to impulse in phase noise')
legend('ct', 'dt')

figure, hold on
impulse(Gs/(Gs+1))
impulse(Gz/(Gz+1))
title('Phase error response to impulse in AWGN noise')
legend('ct', 'dt')

% 
s = tf('s');
z = tf('z', 1/sim.Rs);
figure, hold on
step(1/(s*(Gs+1)))
impulse(z/((z-1)^2*(Gz+1)))
title('Phase error response to ramp in phase noise')
% legend('ct', 'dt')

%
