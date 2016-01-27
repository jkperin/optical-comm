function [Vin,Xout] = QAM_SC_Tx(b,Shifts,Tx,Rs,M)
    
    % Optical Modulator
    Mod = Tx.mod;
    
    % Constellation
    modObj = modem.qammod('M', Tx.const.M, 'SymbolOrder', 'Gray','InputType', 'Bit');
    
    % Sampling rate
    oversamp = 6;     
    Ts       = 1/Rs; 
    NtotTx   = oversamp*length(b);                                           % total number of samples
    tTx      = Ts/oversamp*(0:NtotTx-1);                                     % time vector
    omegaTx  = 2*pi/(Ts/oversamp)/NtotTx*[0:(NtotTx/2-1), (-NtotTx/2):(-1)]; % freq. vector

    Vin = zeros(4,length(tTx));
    Xout = zeros(2*M,length(b));
    B = [];
    for sh = Shifts
        B = [B ; circshift(b,[0,sh])];
    end 
    sum= zeros(1,49146);
    
    for ch = 1:M
        index = 2*log2(Tx.const.M)*(ch-1)+1 : 2*log2(Tx.const.M)*ch;
        X_x = modulate(modObj, B(index(1:length(index)/2),:));
        X_y = modulate(modObj, B(index(length(index)/2+1:end),:));
        
        X = conj([X_x;X_y]);
        Xout((2*ch-1:2*ch),:) = X;
              
        % extract I and Q in x and y
        X1x = real(X(1,:));
        X2x = imag(X(1,:));
        X1y = real(X(2,:));
        X2y = imag(X(2,:));

        % create rectangular waveforms for modulators
        V1x = Tx.Vpi*reshape(repmat(X1x,oversamp,1),oversamp*length(X1x),1)';
        V2x = Tx.Vpi*reshape(repmat(X2x,oversamp,1),oversamp*length(X2x),1)';
        V1y = Tx.Vpi*reshape(repmat(X1y,oversamp,1),oversamp*length(X1y),1)';
        V2y = Tx.Vpi*reshape(repmat(X2y,oversamp,1),oversamp*length(X2y),1)';
       
        % compute composite response of Tx (mod. driver and mod.)
        if strcmp(Tx.Type,'Bessel')
            besfact = [1.000, 1.272, 1.405, 1.515, 1.621];                % factors for cutoff frequencies of Bessel lowpass filters
            [num,denom] = besself(Tx.Ord,2*pi*besfact(Tx.Ord)*Tx.BW);  
            Htx = polyval(num, 1i*omegaTx)./polyval(denom, 1i*omegaTx);        % Bessel LPF
        elseif strcmp(Tx.Type,'Butter');
            [num,denom] = butter(Tx.Ord,2*pi*Tx.BW,'s');
            Htx = polyval(num, 1i*omegaTx)./polyval(denom, 1i*omegaTx);          % Butterworth LPF
        elseif strcmp(Tx.Type,'Cheby');% Chebyshev TypeII 
            R  = 20; 
            Wo = cosh(1/Rx.LPF.Ord*acosh(sqrt(10^(R/10)-1)));
            [num, denom] = cheby2(Rx.LPF.Ord,R,Wo,'s');
            Htx = polyval(num, 1i*omega/(2*pi*Rx.LPF.BW))./polyval(denom, 1i*omega/(2*pi*Rx.LPF.BW));   % case of Chebyshev type II LPF
        end

        % ###################### % 
        deltaomega = omegaTx(2)-omegaTx(1);
        gdtx = -diff(unwrap(phase(Htx)))/deltaomega;                            % group delay at zero frequency
        Htxnd = Htx.*exp(1i*omegaTx*gdtx(1));                                      % freq resp with delay removed
 
        % filter drive waveforms for modulators
        V1x_f = real(ifft(fft(V1x).*Htxnd));
        V2x_f = real(ifft(fft(V2x).*Htxnd)); 
        V1y_f = real(ifft(fft(V1y).*Htxnd)); 
        V2y_f = real(ifft(fft(V2y).*Htxnd));

        Vin_SC = [V1x;V2x;V1y;V2y];
         
        fI = (ch-1)*Rs - (M-1)*Rs/2;
        Vin = Vin + freqshift(Vin_SC,tTx,fI);
        
    end
    
    Vix = real(Vin(1,:)+1i*Vin(2,:));
    Vqx = imag(Vin(1,:)+1i*Vin(2,:));
    Viy = real(Vin(3,:)+1i*Vin(4,:));
    Vqy = imag(Vin(3,:)+1i*Vin(4,:));

    Vix_pd = 2*Tx.Vpi/pi*asin(Vix/max(abs(Vix)));
    Vqx_pd = 2*Tx.Vpi/pi*asin(Vqx/max(abs(Vqx)));
    Viy_pd = 2*Tx.Vpi/pi*asin(Viy/max(abs(Viy)));
    Vqy_pd = 2*Tx.Vpi/pi*asin(Vqy/max(abs(Vqy)));
    
   
    % Frequency response compensation         
    a   = .005*100*sqrt(abs(omegaTx)/2/pi*1e-9);
    Hel = (1-exp(Mod.L*(-a + 1j*Mod.d_12*omegaTx)))./(a-1j*Mod.d_12*omegaTx)/Mod.L;
    Hel(1) = 1;
    Vix = real(ifft(fft(Vix_pd).*(1./Hel)));
    Vqx = real(ifft(fft(Vqx_pd).*(1./Hel))); 
    Viy = real(ifft(fft(Viy_pd).*(1./Hel))); 
    Vqy = real(ifft(fft(Vqy_pd).*(1./Hel))); 
    
    
    ratio = 24/oversamp;
    Vin_ix = interp(Vix,ratio);
    Vin_qx = interp(Vqx,ratio);
    Vin_iy = interp(Viy,ratio);
    Vin_qy = interp(Vqy,ratio);

    Vin = [Vin_ix;Vin_qx;Vin_iy;Vin_qy];
            
end