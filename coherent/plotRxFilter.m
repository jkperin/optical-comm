omega = 2*pi/(Ts/Nstep)/Ntot*(0:Ntot);
% Bessel
besfact = [1.000, 1.272, 1.405, 1.515, 1.621];
[num,denom] = besself(Rx.LPF.Ord,besfact(Rx.LPF.Ord));
Hrx1 = polyval(num, 1i*omega/(2*pi*Rx.LPF.BW))./polyval(denom, 1i*omega/(2*pi*Rx.LPF.BW));   % case of Bessel LPF

% Butterworth
[num,denom] = butter(Rx.LPF.Ord,1,'s');
Hrx2 = polyval(num, 1i*omega/(2*pi*Rx.LPF.BW))./polyval(denom, 1i*omega/(2*pi*Rx.LPF.BW));   % case of Butterworth LPF

% Chebyshev TypeII 
R  = 20; 
Wo = cosh(1/Rx.LPF.Ord*acosh(sqrt(10^(R/10)-1)))
[num, denom] = cheby2(Rx.LPF.Ord,R,Wo,'s');
Hrx3 = polyval(num, 1i*omega/(2*pi*Rx.LPF.BW))./polyval(denom, 1i*omega/(2*pi*Rx.LPF.BW));   % case of Chebyshev type II LPF



plot(omega/(2*pi*1e9),20*log10(abs(Hrx1)),omega/(2*pi*1e9),20*log10(abs(Hrx2))...
    ,omega/(2*pi*1e9),20*log10(abs(Hrx3)),'LineWidth',1.5)
legend('Bessel','Butterworth','Chebyshev TypeII')
xlabel('Frequency (GHz)'); ylabel('Magnitude Response')
title(sprintf('Receiver filter, BW = %.1fR_{s}',Rx.LPF.BW/Rs))
axis([0 80 -30 5])
