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

idx = (size(Hi, 1)+1)/2;
Hic = Hi;
Hqc = Hq;
hic = Hi(idx, :);
hqc = Hq(idx, :);
Hic(idx, :) = [];
Hqc(idx, :) = [];

%% Perform optimization
warning('off', 'MATLAB:nargchk:deprecated')
cvx_begin
    variable pI(Ntaps)
    variable pQ(Ntaps)
%     minimize ( norm(Hic*p, 2) + norm(Hqc*p, 2) - hic*p - hqc*p)
%     minimize ( norm(Hi*p -e, 2) + norm(Hq*p - e, 2) )
%     minimize ( norm(Hic*pI - Hqc*pq, 2) + norm(Hqc*pI + Hic*pq, 2) )
    minimize ( norm(Hi*pI - Hq*pQ - e, 2)  +  norm(Hi*pQ + Hq*pI - e, 2) )
%       minimize (norm(pI + 1j*pQ, 2))
    subject to
          pQ == zeros(Ntaps, 1)
%         pI == [0;  pQ(1:end-1)]
%         pQ == [0; 0; pI(1:end-2)]
%           pI == [pQ(1:end-1); 0]
%         pI == pQ
%           Hic*pI == Hqc*pQ
%           Hic*pQ == -Hqc*pI
cvx_end

if not(strcmpi(cvx_status, 'Solved'))
    warning('Filter optimization failed')
end

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
    stem(n, real((Hi + 1j*Hq)*p))
    stem(n, imag((Hi + 1j*Hq)*p))
%     stem(n, Hi*e, '--')
%     stem(n, Hq*e, '--')
    title('Filtering')
    
    subplot(222), box on, hold on
    stem(-(Ntaps-1)/2:(Ntaps-1)/2, real(p))
    stem(-(Ntaps-1)/2:(Ntaps-1)/2, imag(p))
    title('Filter taps')
    
    P = freqz(p, 1, f, Rs);
    subplot(223)
    plot(f/1e9, abs(P).^2)
    a = axis;
    axis([-Rs/2e9 Rs/2e9 a(3:4)])
    xlabel('Frequency (GHz)')
    ylabel('|P(f)|^2')
    title('Amplitude response')
    
    subplot(224)
    plot(f/1e9, unwrap(angle(P)))
    a = axis;
    axis([-Rs/2e9 Rs/2e9 a(3:4)])
    xlabel('Frequency (GHz)')
    ylabel('arg(P(f))')
    title('Phase response') 
end