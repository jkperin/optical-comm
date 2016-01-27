c      = 3.00e8;
a_m    = 0;
L      = .06;  %cm
n_r    = 2.15; 
hold on
% for ratio  = [.9 .95 .99];
%     a_m    = .0461*100*sqrt(abs(w)/2/pi*1e-9);
%     n_m    = n_r*ratio; 
%     d_12   = (n_m -n_r)/c;
%     w      = (0:.1:600)*1e9;
%     Phi    = (1-exp(-a_m*L + 1j*d_12*w*L))./(a_m-1j*d_12*w);
%     plot(w/2/pi*1e-9,10*log10(abs(Phi)/max(abs(Phi))),'b--');
%     axis([0 50 -10 1])
% end
for L = [0.05]
    ratio = .95;
    n_m    = n_r*ratio; 
    d_12   = (n_m -n_r)/c;
    w      = (-600:.1:600)*1e9;
    a_m    = .01*100*sqrt(abs(w)/2/pi*1e-9);
    Phi    = (1-exp(-a_m*L + 1j*d_12*w*L))./(a_m-1j*d_12*w);
    plot(w/2/pi*1e-9,10*log10(abs(Phi)/max(abs(Phi))),'r-');
    %axis([0 50 -10 1])
end
plot([0,50],[-3 -3],'k--')
hold off
%axis([0 50 -5 0])