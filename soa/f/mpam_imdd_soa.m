function ber = mpam_imdd_soa(mpam, tx, soa, rx, sim)

% Generate data
dataTX = randi([0 mpam.M-1], [sim.Nsymb 1]);

x = mpam.a(gray2bin(dataTX, 'pam', mpam.M) + 1);

% Transmitter pulse shape
gtx = ones(1, sim.Mct);

% Continuous time
xt = kron(x, gtx);

xt = reshape(xt.', 1, sim.Nsymb*sim.Mct);

% Scale so that average power is equal to the received power
P = xt;

%% RIN noise
% nrin = 0;
% if sim.rin
%     % Calculate double-sided RIN PSD, which depends on the average power
%     Srin = 10^(tx.RIN/10)*(rx.Gapd*rx.R*P).^2;
%     
%     Bw = sim.Mct*mpam.Rb/log2(mpam.M);
% 
%     nrin = sqrt(Srin*2*Bw).*randn(size(P));
% end

%% Electric signal without chirp
Et = sqrt(P);

%% Add noise due to SOA
if sim.soa  
    Bw = sim.Mct*mpam.Rs;
    
    nr =sqrt(soa.Seq(soa.G)*Bw/2)*randn(size(Et));
    ni =sqrt(soa.Seq(soa.G)*Bw/2)*randn(size(Et));
    
    Et = sqrt(soa.G)*Et + nr + 1j*ni;
end

%% Bandpass filter
% dt = 1/(sim.Mct*mpam.Rs);
% df = 1/(dt*length(Et));
% fs = 1/dt;
% f = -fs/2:df:fs/2-df;
% imp_length = 256;
% 
% [bfilt, afilt] = design_filter('gaussian', 4, mpam.Rs/fs, 1);
% 
% gbp = impz(bfilt, afilt, imp_length).';    
% 
% gbp = gbp/sum(gbp);
% gbp_delay = grpdelay(bfilt, afilt, 1);    % calculates the group delay by approximating the filter as an FIR filter
%                                             % whose impulse response is given by tx.gdac   
% 
% % Remove group delay
% Gbp = fft(gbp, length(f));
% Gbp = fftshift(Gbp).*exp(1j*2*pi*f*gbp_delay/fs); 
% 
% figure(1), hold on
% Ef = fftshift(abs(fft(Et-mean(Et))));
% plot(f, Ef/max(Ef))
% plot(f, abs(Gbp), 'k')

% % Filter optical signal before direct detection
% Et = ifft(fft(Et).*ifftshift(Gbp));

% Square law detection
Irx = rx.R*abs(Et).^2;

% % Matched filter and Donwsample
% Irx = filter(fliplr(gtx), 1, Irx);
% 
% Irx = Irx(sim.Mct:sim.Mct:end);

% Add thermal noise from TIA
Idet = Irx + sqrt(rx.N0*mpam.Rs)*randn(size(Irx));

% Rescale signal for detection
y = tx.Prec/mean(Idet)*Idet.';

% Demodulate
dataRX = sum(bsxfun(@ge, y, mpam.b.'), 2);
dataRX = bin2gray(dataRX, 'pam', mpam.M);

% True BER
Ndisc = 20;
[~, ber] = biterr(dataRX(Ndisc+1:end-Ndisc), dataTX(Ndisc+1:end-Ndisc));

% BER using approximation BER = SER/log2(M)
% ber = mean(dataTX ~= dataRX)/log2(mpam.M);
