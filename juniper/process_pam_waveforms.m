%% Process scope data
function ber_count = process_pam_waveforms(file, dacfile)

addpath f/
addpath ../f/
addpath ../mpam/

% file = 'data/waveforms/b2b_OSNR/wave-4PAM-56Gbps-preemph-0km-osnr=42dB.h5';
% dacfile = 'data/waveforms/b2b_OSNR/pam4_rect_Rb=55Gbps_preemph';
Dac = load(dacfile);

%
fs = 80e9; % Receiver sampling rate
ros = Dac.sim.ros.txDSP; % DAC oversampling rate
mpam = Dac.mpam;

%% Load channels
hinfo =  hdf5info(file);
ch1 = h5read(file, '/Waveforms/Channel 1/Channel 1Data');
ch2 = h5read(file, '/Waveforms/Channel 2/Channel 2Data');
ch3 = h5read(file, '/Waveforms/Channel 3/Channel 3Data');
ch4 = h5read(file, '/Waveforms/Channel 4/Channel 4Data');

ch1 = double(ch1).';
ch2 = double(ch2).';
ch3 = double(ch3).';
ch4 = double(ch4).';

%% Orthornomalization
ch1 = ch1 - mean(ch1);
ch1 = ch1/sqrt(mean(abs(ch1).^2));
ch2 = ch2 - mean(ch2);
ch2 = ch2/sqrt(mean(abs(ch2).^2));
ch3 = ch3 - mean(ch3);
ch3 = ch3/sqrt(mean(abs(ch3).^2));
ch4 = ch4 - mean(ch4);
ch4 = ch4/sqrt(mean(abs(ch4).^2));

ch2 = ch2 - mean(ch1.*ch2)*ch1;
ch4 = ch4 - mean(ch3.*ch4)*ch4;

%% Pre-filtering
f = freq_time(length(ch1), fs);
Filt = design_filter('butter', 5, mpam.Rs/(fs/2)); % half of the sampling rate

ch1 = real(ifft(fft(ch1).*ifftshift(Filt.H(f/fs))));
ch2 = real(ifft(fft(ch2).*ifftshift(Filt.H(f/fs))));
ch3 = real(ifft(fft(ch3).*ifftshift(Filt.H(f/fs))));
ch4 = real(ifft(fft(ch4).*ifftshift(Filt.H(f/fs))));

X = ch1 + 1j*ch2;
Y = ch3 + 1j*ch4;

%% Frequency offset compensaion
% [~, foff_ind] = max(abs(fftshift(fft(X))));
% foff = f(foff_ind)
% 
% X = X.*exp(-1j*2*pi*foff*t);
% Y = Y.*exp(-1j*2*pi*foff*t);
% 
% figure, plot(f, abs(fftshift(fft(X).^2)))
% figure, plot(t, real(X), t, imag(X))

%% Convert to intensity
Xrx = abs(X).^2 + abs(Y).^2; % Get intensity in X-pol

%% Antialiasing filter and reduce noise
f = freq_time(length(Xrx), fs);
Filt = design_filter('butter', 5, 0.7*mpam.Rs/(fs/2)); % half of the sampling rate
Xfilt = real(ifft(fft(Xrx).*ifftshift(Filt.H(f/fs))));

%% Resample
% Resample in order to have integer number of samples per symbol specified by ros
[p, q] = rat(ros*mpam.Rs/fs);
Xrx = resample(Xfilt, p, q); % resample to match symbol rate at the DAC

Xrx = Xrx - mean(Xrx);

xref = mpam.signal(Dac.dataTX); % reference signal
xref = xref - mean(xref);

% Align received signal
[c, lags] = xcorr(Xrx, xref);
cmax = max(abs(c));
threshold = 0.7;
ind = find(abs(c) >= cmax*threshold & lags > 1);
pos = lags(ind);
remapping =(sign(c(ind(1))) == -1);

figure, plot(lags, c)

Xrx = Xrx(pos(1):end);
Npatterns = floor(length(Xrx)/length(xref));
N = Npatterns*length(xref);
Xrx = Xrx(1:N);

%% PAM symbol detection
yk = Xrx;

%% Finer gain control
yk = yk - mean(yk);
yk = yk*sqrt(var(mpam.a)/(mean(abs(yk).^2)));
% yk = yk + mean(mpam.a);

% f = freq_time(length(yk), mpam.Rs*ros)/(mpam.Rs*ros);
% yk = real(ifft(fft(yk).*ifftshift(exp(1j*2*pi*f*1/2))));

mpam.b = mpam.b - mean(mpam.a);
mpam.a = mpam.a - mean(mpam.a);

%% Equalization
% Resample to have at least 2 samples per symbol
[p, q] = rat(2/ros);
yk = resample(yk, p, q); 

eq.type = 'Adaptive TD-LE';
eq.Ntaps = 15;
eq.mu = 1e-2;
eq.ros = 2;
eq.Ntrain = Inf;
eq.Ndiscard = [1e5 1024];

if remapping 
    remap = [2 3 0 1];
    dataTX = remap(Dac.dataTX + 1);
else
    dataTX = Dac.dataTX;
end
dataTX = repmat(dataTX, 1, Npatterns);
eq.trainSeq = dataTX;

% xref2 = real(ifft(fft(2*mpam.signal(dataTX)).*ifftshift(exp(1j*2*pi*f*9/2))));
% xref3 = repmat(Dac.xt(1:Dac.sim.Mct/2:end), 1, Npatterns);
% xref3 = xref3 - mean(xref3);
% xref3 = xref3*sqrt(var(mpam.a)/(mean(abs(xref3).^2)));
% xref3 = xref3 + mean(mpam.a);
sim.Nsymb = length(eq.trainSeq);
sim.Ndiscard = 512;
sim.f = freq_time(length(yk), mpam.Rs*eq.ros);
[yd, eq] = equalize(eq, yk, [], mpam, sim, true);

% Symbols to be discard in BER calculation
ndiscard = [1:eq.Ndiscard(1)+sim.Ndiscard (sim.Nsymb-eq.Ndiscard(2)-sim.Ndiscard):sim.Nsymb];
ydfull = yd;
yd(ndiscard) = []; 
dataTX(ndiscard) = [];

%% Demodulate
dataRX = mpam.demod(yd.');

%% Counted BER
[~, ber_count] = biterr(dataRX, dataTX)

% mpam = mpam.norm_levels();
figure(101), clf, box on, hold on
h1 = plot(ydfull, 'o');
a = axis;
h2= plot(a(1:2), (mpam.a*[1 1]).', '-k');
h3 = plot(a(1:2), (mpam.b*[1 1]).', '--k');
h4 = plot(eq.Ndiscard(1)*[1 1], a(3:4), ':k');
h5 = plot((sim.Nsymb-eq.Ndiscard(2))*[1 1], a(3:4), ':k');
legend([h1 h2(1) h3(1) h4], {'Equalized samples', 'PAM levels',...
    'Decision thresholds', 'BER measurement window'})
title('Symbol-rate samples after equalization')
% axis([1 sim.Nsymb -0.2 1.2])
drawnow

%% 4-PAM Bias analysis
if mpam.M == 4
	p(1) = mean(yd(yd < mpam.b(1)));
    p(2) = mean(yd(yd > mpam.b(1) & yd < mpam.b(2)));
    p(3) = mean(yd(yd > mpam.b(2) & yd < mpam.b(3)));
    p(4) = mean(yd(yd > mpam.b(3)));
    
    mzm_nonlinearity = @(levels, V) abs(sin(pi*(levels*V(1) + V(2)))).^2;
    
    [Vset, fval, exitflag] = fminsearch(@(V) norm(p.' - (mzm_nonlinearity(mpam.a, V) - mean(mzm_nonlinearity(mpam.a, V)))), [0.5 0.5]);
    
    if exitflag ~= 1
        disp('4-PAM bias control did not converge')
    end
    
    fprintf('Vgain/Vpi = %.2f\n', Vset(1))
    fprintf('Vbias/Vpi = %.2f\n', Vset(2))
    
    figure(102), clf, box on, hold on
    h1 = plot([1 2], ((mpam.a + 0.5)*[1 1]), '-k');
    h2 = plot([1 2], (mzm_nonlinearity(mpam.b, Vset)*[1 1]).', '--k');
    h3 = plot([1 2], (mzm_nonlinearity(mpam.a, Vset)*[1 1]), ':k');
    legend([h1(1) h2(1) h3(1)], {'Desired PAM levels', 'Decision thresholds', 'Distorted PAM levels'})
    title('4-PAM Bias analysis')
    % axis([1 sim.Nsymb -0.2 1.2])
    drawnow
end
    


    
    
    



% S = comm.SymbolSynchronizer('TimingErrorDetector', 'Gardner (non-data-aided)', 'SamplesPerSymbol', 4);
% yd = step(S, yk.');
