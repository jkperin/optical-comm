%% Compare performance of M-PAM in 3 different scenarios
% 1. PIN receiver, no amplifier
% 2. SOA & PIN receiver
% 3. APD
clear, clc, close all

%% Simulation parameters
sim.Nsymb = 2^14; % Number of symbols in montecarlo simulation
sim.Mct = 16;     % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

sim.ros = 2;

%% M-PAM
mpam = PAM(4, 100e9, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

% Modulator frequency response
tx.modulator.fc = 30e9; % modulator cut off frequency
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds
tx.modulator.H = @(f) (1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2)).*exp(1j*2*pi*f*tx.modulator.grpdelay);

%% Fiber
fiber = fiber(); % fiber(L, att(lamb), D(lamb))

%% Receiver
rx.N0 = (20e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
rx.antialiasing = design_filter('bessel', 4, mpam.Rs/(sim.fs/2));

dataTX = randi([0 mpam.M-1], 1, sim.Nsymb);
x = mpam.mod(dataTX, sim.Mct);

p = real(ifft(fft(x).*ifftshift(tx.modulator.H(f))));

p = p + rx.N0*sim.fs/2*randn(size(p));

pa = real(ifft(fft(p).*ifftshift(rx.antialiasing.H(f/sim.fs))));

figure, hold on
plot(t, x)
plot(t, pa)

pk = pa(1:sim.Mct/sim.ros:end);
tk = t(1:sim.Mct/sim.ros:end);

% pk = resample(pa, sim.ros, sim.Mct);
% tk = resample(sim.t, sim.ros, sim.Mct);

plot(tk, pk, 'o')

% Time-domain equalizer
% pk = p((Mct-1)/2:Mct:end);
Ntaps = 17;
Ntrain = 5e3;
W = [zeros(Ntaps-1, 1); 1];
y = zeros(size(pk));
b = mpam.mod(dataTX, 1);
mu = 1e-2;
n = 1;
for k = Ntaps:length(pk)
    z = pk(k-Ntaps+1:k);
    
    y(k) = sum(W.*z);
    
    if mod(k, sim.ros) == 0
        if n < Ntrain % Training
%             e((k-1)/sim.ros) = y(k) - b((k-1)/sim.ros);

            e(k/sim.ros) = y(k) - b(k/sim.ros-1);
        
        else
%             e((k-1)/sim.ros) = y(k) - mpam.mod(mpam.demod(y(k)), 1);
            
            e(k/sim.ros) = y(k) - mpam.mod(mpam.demod(y(k)), 1);
            
        end
        n = n + 1;
        
%         W = W - 2*mu*e((k-1)/sim.ros)*z; 
        W = W - 2*mu*e(k/sim.ros)*z;
    end
end

% plot(tk, b, 'o')

plot(tk(2:sim.ros:end), y(2:sim.ros:end), 'o')
        

figure, plot(e)
        
        
H = freqz(W, 1, sim.f, sim.fs);
figure
plot(sim.f/1e9, abs(H).^2)


%% Equalizer
% Nfft = sim.Nsymb*sim.ros;
% ff = -0.5:1/Nfft:0.5-1/Nfft;
% ef = design_filter('matched', mpam.pshape, 1/sim.ros);
% % W = 1./(tx.modulator.H(ff*mpam.Rs*sim.ros).*ef.H(ff));
% W = 1;
% % figure
% % plot(
% 
% Pk = fftshift(fft(pk));
% Pk = W.'.*Pk;
% 
% Ns = round((sim.ros-1)/sim.ros*length(pk))/2;
% % Pk = [0:length(pk)/2 length(pk)/2-1:-1:1];
% 
% Pk2 = Pk(Ns+1:end-Ns);
% Pk2(1:Ns+1) = Pk2(1:Ns+1) + conj(fliplr(Pk(1:Ns+1)));
% Pk2(end-Ns+2:end) = Pk2(end-Ns+2:end) + conj(fliplr(Pk(end-Ns+2:end)));
% 
% figure, hold on
% plot(-128:127, abs(Pk))
% plot(-64:63, abs(Pk2))
% 
% yk = real(ifft(ifftshift(Pk2)));
% 
% figure(1)
% plot(tk(1:sim.ros:end), yk/2, '-*')


