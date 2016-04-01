function W = Weiner()
    close all
    l0 = 12;
    k = 50;
    N0 = 2.0768e-016;
    [a_ii ~] = Ql(l0,l0,k,k,0);
    a_ii = conj(a_ii');
    a = [a_ii,zeros(size(a_ii));zeros(size(a_ii)),a_ii];
    A_ii = zeros(size(a_ii,1),size(a_ii,1));
    for l=1:24;
        for n=(k-5):(k+5)
            [q,p] = Ql(l0,l,k,n,0);
            A_ii = A_ii + conj(q')*q + N0*corr(p,size(a_ii,1));
        end
    end
    A = [A_ii,zeros(size(A_ii));zeros(size(A_ii)),A_ii];
    W = inv(A)*a;
end

function [Samples,P] = Ql(l0,l,k,n,plots)
    global Ts Nstep oversamp;
    % ============ CONSTANTS ============
    lambda = 1550e-9;                                                           % wavelength (m)
    c      = 3.00e8;                                                            % speed of light (m/s)                                                        
    % ======== Co-OFDM parameteres =======
    Rs     = 12.5e9;                                                            % Symbol rate per carrier (symbol/s)
    Ts     = 1/Rs;                                                              % Symbol interval per carrier (s)

    % ====== Simulation parameter ======
    % For modulation, propagation and detection of carriers -------------------
    Nsymb = 100;
    Nstep = 96;                                                                 % Time steps per symbol interval (mulitple of oversamp)
    Ntot  = Nstep*Nsymb;                                                        % total number of time steps in waveforms                
    t     = (0:Ntot-1)*Ts/Nstep;                                                % Continuous time vector
    omega = 2*pi/(Ts/Nstep)/Ntot*[0:(Ntot/2-1), (-Ntot/2):(-1)];                % continuous freq. vector
    % ====================================
    b0 = (t>=(k-1)*Ts) .* (t<k*Ts); B0 = abs(fft(b0));
    % ----- b(t) ----- %
    Tx.Ord = 5;
    Tx.BW    = 1.2*Rs;
    besfact = [1.000, 1.272, 1.405, 1.515, 1.621];                              % factors for cutoff frequencies of Bessel lowpass filters
    [num,denom] = besself(Tx.Ord,2*pi*besfact(Tx.Ord)*Tx.BW);  
    Htx = polyval(num, 1i*omega)./polyval(denom, 1i*omega);                     % Bessel LPF
    deltaomega = omega(2)-omega(1);
    gdtx = -diff(unwrap(phase(Htx)))/deltaomega;                                % group delay at zero frequency
    Htxnd = Htx.*exp(1i*omega*gdtx(1));                                         % freq resp with delay removed                            
    b = real(ifft(fft(b0).*Htxnd));  B = abs(fft(b));
    % -------- bl(t) ---------- %
    %l0 = 12;
    %l  = 12;
    df_l = (l-l0)*Rs;
    b_l = freqshift(b,t,df_l);    Bl = abs(fft(b_l));
    % ====== Fiber (excluding PMD) ========
    Fiber.Dpsnmkm = 17;                                                         % SMF CD (ps/(nm*km))
    Dsmm          = 1e-6*Fiber.Dpsnmkm;           
    Fiber.beta2   = -lambda^2/(2*pi*c)*Dsmm;                                    % SMF beta2
    Fiber.alphadBkm = 0-0.25;                                                    % SMF loss (dBm/km)
    Fiber.Lkm     = 80;                                                         % SMF span length (km)
    Fiber.Tf      = 10^(Fiber.alphadBkm*Fiber.Lkm/10);                          % Power transmittance of fiber segment
    Lseg          = 1e3*Fiber.Lkm;                                              % Fiber segment length (m)
    Fiber.H       = sqrt(Fiber.Tf)*exp(-1i*Fiber.beta2/2*Lseg*omega.^2);        % Freq resp of fiber segment (excluding PMD)
    c = ifft(fft(b_l).*Fiber.H);  C = abs(fft(c));

    % ====== Anti-aliasing filter ====== %
    Rx.LPF.BW   = 2.4*Rs;                                                       % Bandwidth at -3 dB of rec. filt.
    Rx.LPF.Ord  = 5;   
    [num,denom] = butter(Rx.LPF.Ord,2*pi*Rx.LPF.BW,'s');
    Hrx = polyval(num, 1i*omega)./polyval(denom, 1i*omega);  
    deltaomega = omega(2)-omega(1);
    gdrx = -diff(unwrap(phase(Hrx)))/deltaomega;                               % group delay at zero frequency
    Hrxnd = Hrx.*exp(1i*omega*gdrx(1));                                         % freq resp with delay removed
    q1 = (ifft(fft(c).*Hrxnd));  Q1 = abs(fft(q1));
    % ------- CD Compensation -------- %                                                     % Window type. Options include: bartlett, chebwin, hamming, hanning, rectwin (see MATLAB "help window"). 
    Ltot    = 1e3*Fiber.Lkm;                                                    % total SMF length (m)
    Hcdcomp = exp(1i*Fiber.beta2/2*Ltot*omega.^2);  
    q2 = ifft(fft(q1).*Hcdcomp); Q2 = abs(fft(q2));
    % ------- Carrier Separation ------- % 
    HT2nd = 0.5*(exp(1i*omega*Ts/4)+exp(-1i*omega*Ts/4));                 % delay by T/2 and add (shifted to give zero group delay)
    HT4nd = 0.5*(exp(1i*omega*Ts/8)+exp(-1i*omega*Ts/8));                 % delay by T/4 and add (shifted to give zero group delay)
    Hcarsep = HT2nd.*HT2nd.*HT4nd;                                              % overall carrier separation filter
    q3 = ifft(fft(q2).*Hcarsep); Q3 = abs(fft(q3));
    % ---- Sampler -----% 
    oversamp = 6;
    firstsamp = 1 + floor(0.5*Nstep);                                     % index of first sample
    sampincr = floor(Nstep/oversamp);
    lastsamp = floor(0.5*Nstep) + oversamp*Nsymb*floor(Nstep/oversamp);   % index of first sample
    sampind = mod(firstsamp:sampincr:lastsamp,Nsymb*Nstep);          % indices of samples in x pol. (modulo wraps around)
    qs = q3(sampind);           % samples of inphase in x pol.
    tsamp = t(sampind); 

    Lfilt = 10;
    eqindk = mod(floor(oversamp*(n-1))+Lfilt:-1:floor(oversamp*(n-1))-Lfilt,length(tsamp))+1;
    Samples = qs(eqindk);

    % ---- Noise term ---- %
    P = fftshift(ifft(Hrxnd.*Hcdcomp.*Hcarsep));

    if (plots)
    % =========================================
    figure()
    plot(t,P,'r');
    figure()
    subplot(121)
    plot(t,b0,'r',t,b,'b',t,b_l,'k','LineWidth',1.5);
    axis([(k-1.5)*Ts,(k+.5)*Ts, -1.5 1.5])
    subplot(122)
    plot(omega/2/pi/1e9,10*log10(B0),'r',...
        omega/2/pi/1e9,10*log10(B),'b',...
        omega/2/pi/1e9,10*log10(Bl),'k','LineWidth',1.5)
    axis([-5*Rs/1e9 5*Rs/1e9 10 30]);
    % =========================================
    figure()
    subplot(121)
    plot(t,b_l,'b',t,c,'g',t,q1,'r--',t,q2,'m--','LineWidth',1.5);
    axis([(k-5)*Ts,(k+5)*Ts, -1.5 1.5])
    subplot(122)
    plot(omega/2/pi/1e9,10*log10(Bl),'b',...
        omega/2/pi/1e9,10*log10(C),'g--',...
        omega/2/pi/1e9,10*log10(Hrxnd),'k--',...
        omega/2/pi/1e9,10*log10(Q1),'r',...
        omega/2/pi/1e9,10*log10(Q2),'m--','LineWidth',1.5);
    axis([-5*Rs/1e9 5*Rs/1e9 -5 20]);
    %==========================================
    figure()
    subplot(121)
    plot(t,b_l,'b',t,q3,'m','LineWidth',1.5);
    axis([(k-2.5)*Ts,(k+1.5)*Ts, -1.5 1.5])
    subplot(122)
    plot(omega/2/pi/1e9,10*log10(Q2),'b',...
        omega/2/pi/1e9,10*log10(10*Hcarsep),'k.',...
        omega/2/pi/1e9,10*log10(Q3),'m','MarkerSize',5,'LineWidth',1.5);
    axis([-5*Rs/1e9 5*Rs/1e9 -5 20]);
    %===========================================
    figure()
    plot(t/Ts,q3,tsamp(eqindk)/Ts,qs(eqindk),'ro','LineWidth',1.5)
    axis([min(tsamp(eqindk))/Ts-1,max(tsamp(eqindk))/Ts+1, -1.5 1.5])
    end
end
function c = corr(x,N)
    global Ts Nstep oversamp;
    c = zeros(N,N);
    for i=N
        for j=N     
            c(i,j) = Ts/Nstep*sum(x.*circshift(x,(i-j)*Nstep/oversamp));
        end
    end
end
