function p = design_predistortion_filter(Ntaps, ros, mpam, Fiber, wavelength, verbose)
if mod(Ntaps, 2) == 0
    Ntaps = Ntaps + 1; % make it odd
end

% 
Rs = mpam.Rs;
Mct = ros;
N = 2^12;
fs = Mct*Rs;
f = freq_time(N, fs);

% HtxPshape = sinc(f/Rs);
H = Fiber.Hdisp(f, wavelength);
% H = Fiber.Hdisp(f, wavelength) ;
h = fftshift(ifft(ifftshift(H)));
hi = real(h).';
hq = imag(h).';

n = (-N/2:N/2-1).';
% idx = find(mod(abs(n), Mct) == 0 & abs(n) <= Ntaps*Mct); 
idx = find(abs(n) <= Ntaps*Mct); 
nd = n(idx);
hid = hi(idx); 
hqd = hq(idx); 

figure(343), hold on
hline(1) = plot(n, hi);
hline(2) = plot(n, hq);
stem(nd, hid, 'Color', get(hline(1), 'Color'));
stem(nd, hqd, 'Color', get(hline(2), 'Color'));
plot(n, 1/fs*real(Fiber.hdisp(n/fs, wavelength)), '--')
plot(n, 1/fs*imag(Fiber.hdisp(n/fs, wavelength)), '--')
a = axis;
axis([nd(1) nd(end) a(3:4)])

Hi = toeplitz([hid; zeros(Ntaps-1, 1)], [hid(1) zeros(1, Ntaps-1)]);
Hq = toeplitz([hqd; zeros(Ntaps-1, 1)], [hqd(1) zeros(1, Ntaps-1)]);

% ZF condition
e = zeros(size(Hi, 1), 1); 
e((size(Hi, 1)+1)/2) = 1;
% e((size(Hi, 1)+1)/2-1) = 1;
% e((size(Hi, 1)+1)/2-2) = 1;
% e((size(Hi, 1)+1)/2-3) = 1;

idx = (size(Hi, 1)+1)/2;
Hic = Hi;
Hqc = Hq;
hic = Hi(idx, :);
hqc = Hq(idx, :);
Hic(idx, :) = [];
Hqc(idx, :) = [];

%% Perform optimization
warning('off', 'MATLAB:nargchk:deprecated')
%% I/Q impulse
cvx_begin quiet
    variable pI(Ntaps)
    variable pQ(Ntaps)
    minimize ( norm(Hi*pI - Hq*pQ - e, 2)  +  norm(Hi*pQ + Hq*pI - e, 2) )
    subject to
          pQ == zeros(Ntaps, 1)
%         pI == [0;  pQ(1:end-1)]
%         pQ == [0; 0; pI(1:end-2)]
%           pI == [pQ(1:end-1); 0]
%         pI == -0.1*pQ
%           Hic*pI == Hqc*pQ
%           Hic*pQ == -Hqc*pI
cvx_end

if not(strcmpi(cvx_status, 'Solved'))
    warning('Filter optimization failed')
end

prob(1) = cvx_optval;
probpI{1} = pI;
probpQ{1} = pQ;

%% I impulse, Q zeros
cvx_begin quiet
    variable pI(Ntaps)
    variable pQ(Ntaps)
    minimize ( norm(Hi*pI - Hq*pQ - e, 2)  +  norm(Hi*pQ + Hq*pI, 2) )
    subject to
%           pQ == zeros(Ntaps, 1)
%         pI == [0;  pQ(1:end-1)]
%         pQ == [0; 0; pI(1:end-2)]
%           pI == [pQ(1:end-1); 0]
%         pI == -0.1*pQ
%           Hic*pI == Hqc*pQ
%           Hic*pQ == -Hqc*pI
cvx_end

if not(strcmpi(cvx_status, 'Solved'))
    warning('Filter optimization failed')
end

prob(2) = cvx_optval;
probpI{2} = pI;
probpQ{2} = pQ;

%% Q impulse, I zeros
cvx_begin quiet
    variable pI(Ntaps)
    variable pQ(Ntaps)
    minimize ( norm(Hi*pI - Hq*pQ, 2)  +  norm(Hi*pQ + Hq*pI - e, 2) )
    subject to
%           pQ == zeros(Ntaps, 1)
%         pI == [0;  pQ(1:end-1)]
%         pQ == [0; 0; pI(1:end-2)]
%           pI == [pQ(1:end-1); 0]
%         pI == -0.1*pQ
%           Hic*pI == Hqc*pQ
%           Hic*pQ == -Hqc*pI
cvx_end

if not(strcmpi(cvx_status, 'Solved'))
    warning('Filter optimization failed')
end

prob(3) = cvx_optval;
probpI{3} = pI;
probpQ{3} = pQ;

%% 
prob
[minprob, idx] = min(prob)
pI = probpI{idx};
pQ = probpQ{idx};

%%
pref = hid - 1j*hqd;
pref = pref/abs(sum(pref));

p = pI + 1j*pQ;
% p  = pI;
p = p/abs(sum(p));
% p = pinv(Hi+Hq)*e;
% p = p/abs(sum(p));

if exist('verbose', 'var') && verbose
    figure(454), clf
    n = -(size(Hi, 1)-1)/2:(size(Hi, 1)-1)/2;
    e = zeros(size(Hi, 2));
    e((size(Hi, 2)+1)/2) = 1;
    subplot(221), hold on, box on
    stem(n, real((Hi + 1j*Hq)*p), 'Linewidth', 2)
    stem(n, imag((Hi + 1j*Hq)*p), 'Linewidth', 2)
%     stem(n, Hi*e, '--')
%     stem(n, Hq*e, '--')
    xlabel('Sample', 'FontSize', 18)
    ylabel('Filter output', 'FontSize', 18)
    title('Filtering')
    set(gca, 'FontSize', 18)
    a = axis;
    axis([n(1) n(end) a(3:4)])
    
    subplot(222), box on, hold on
    stem(-(Ntaps-1)/2:(Ntaps-1)/2, real(p), 'Linewidth', 2)
    stem(-(Ntaps-1)/2:(Ntaps-1)/2, imag(p), 'Linewidth', 2)
    xlabel('Sample', 'FontSize', 18)
    ylabel('Filter taps', 'FontSize', 18)
    set(gca, 'FontSize', 18)
    title('Filter taps (real & imag)')
    a = axis;
    axis([-(Ntaps-1)/2 (Ntaps-1)/2 a(3:4)])
    
    P = freqz(p, 1, f, Rs);
    subplot(223)
    plot(f/1e9, abs(P).^2, 'Linewidth', 2)
    a = axis;
    axis([-Rs/2e9 Rs/2e9 a(3:4)])
    xlabel('Frequency (GHz)', 'FontSize', 18)
    ylabel('|P(f)|^2', 'FontSize', 18)
    title('Amplitude response')
    set(gca, 'FontSize', 18)
    
    subplot(224)
    plot(f/1e9, unwrap(angle(P)), 'Linewidth', 2)
    a = axis;
    axis([-Rs/2e9 Rs/2e9 a(3:4)])
    xlabel('Frequency (GHz)', 'FontSize', 18)
    ylabel('arg(P(f))', 'FontSize', 18)
    title('Phase response') 
    set(gca, 'FontSize', 18)
end