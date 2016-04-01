function [xhat, thetahat, Delta] = feedforward_cpr(y, Cpr, sim)
% CPR.m   Milad Sharif  02-22-13
% Simulates carrier recovery for QPSK in two polarizations.
% See: E. Ip and J. M. Kahn, "Feedforward Carrier Recovery for Coherent Optical Communications",
% J. of Lightwave Technol., vol. 25, no. 9, pp. 2675-2692, September 2007. 
% and: E. Ip, A. P. T. Lau, D. J. F. Barros and J. M. Kahn, "Coherent Detection in
% Optical aFiber Systems", Optics Express, vol. 16, no. 2, pp. 753-791,
% January 21, 2008.

% Inputs:
% y : input signal in two polarizations
% Cpr : sctruct containing carrier phase recovery parameters
% Cpr.{phaseEstimation, varPN, SNRdB (optional), Ntaps, FilterType}

    y1 = y(1,:).';
    y2 = y(2,:).';
    
    M = sim.M; % QAM order (must be 4 if Cpr.type ~= 'DD')
    Ntrain = Cpr.Ntrain;
    Ytrain = Cpr.Ytrain;
    sigmapsq = Cpr.varPN;
    crtype   = Cpr.phaseEstimation;    
    L     = Cpr.Ntaps; % FIR: actual length of Whd (for DD or NDA 1) or of Wsd (for NDA 2) employed, IIR: length of Whd (for DD or NDA 1) or of Wsd (for NDA 2) employed in approximating the IIR filter
    filttype = Cpr.FilterType;

    ndtype = '2';                   % type of carrier recovery if nda.
                                    % '1' = as in reference quoted above.
                                    % '2' = reverse ordering of fourth power and digital filtering.

    % perform some computations related to carrier recovery
    probx = 1; % symbols probability (only ~= 1 when shapping is used)
    x = qammod(0:M-1, M); % QAM symbols
    eta_c = mean(abs(x.^2).*probx) * mean((abs(x).^(-2)).*probx);   % constellation penalty

    if isfield(Cpr, 'SNRdB') % SNR was specified
        snr_des_dB = Cpr.SNRdB;
    else
        snr_des_notrack_dB = 8         % SNR per symbol at which the carrier recovery is optimized when snrtrack = 'n'
        snr_des_dB = snr_des_notrack_dB;
    end

    snr_des = 10^(0.1*snr_des_dB);
    
    if strcmp(crtype,'DD')
        sigmandsq = 0.5*eta_c/snr_des;                                          % variance of n'_k
    elseif strcmp(crtype,'NDA')
        sigmandsq = (1/(2*M^2))*sum(((gamma(M+1)./(gamma([1:M]+1).*gamma(M-[1:M]+1))).^2).*gamma([1:M]+1)./(snr_des.^[1:M]));
    end

    if strcmp(filttype,'FIR')
        Delta = floor(0.5*(L-1));                           % delay of Whd (for DD or NDA 1) or of Wsd (for NDA 2)
        if strcmp(crtype,'DD')
            d = floor(0.5*(L-1));                       % length of Wsd for DD
        end
    elseif strcmp(filttype,'IIR')
        Delta = 0;                           % delay of Whd (for DD or NDA 1) or of Wsd (for NDA 2)
        if strcmp(crtype,'DD')
            d = floor(0.5*(L-1));                       % length of Wsd for DD
        end
    end


    N = size(y,2);
    thetahat = zeros(N,2);                                                  % hard-decision phase estimate
    x1hat = zeros(N,1);                                                     % symbols after carrier de-rotation
    x2hat = zeros(N,1);
    x1D = zeros(N,1);                                                       % detected symbols
    x2D = zeros(N,1);
    
    % Computer Wiener filters and perform carrier recovery
    if strcmp(crtype,'DD')
        Wsd = wiener_fir(d,0,sigmandsq,sigmapsq);
        Whd = wiener_fir(L,Delta,sigmandsq,sigmapsq);
       
        psi = zeros(2,N);                                                   % soft-decision phase estimate
        theta_tilde = zeros(2,N);                                           % input to soft-decision phase estimator
        for k=L:N

            % Make soft decision on symbol y(k)
            yk = y(:,k);
            sk = yk.*exp(-1i*theta_tilde(:,k));
            if Ntrain < k % training sequence
                xhatk_sd = Ytrain(:, k);
            else % switch to decision directed
                xhatk_sd = detect(sk);
            end
            psik_tilde = angle(yk)-angle(xhatk_sd);
            p  = floor((psi(:,k-1)-psik_tilde+pi)/(2*pi));
            psi(:,k) = psik_tilde + p*(2*pi);
            theta_tilde(:,k+1) = psi(:,k:-1:k-d+1) * Wsd;
            thetahat(k-Delta, :) = psi(:,k:-1:k-L+1) * conj(Whd);
            % Detect symbol
            x1hat(k-Delta) = y1(k-Delta)*exp(-1i*thetahat(k-Delta, 1));
            x2hat(k-Delta) = y2(k-Delta)*exp(-1i*thetahat(k-Delta, 2));
        end;

    elseif strcmp(crtype,'NDA')
        if ndtype == '1';
            Whd = wiener_fir(L,Delta,sigmandsq,sigmapsq);

            psi = zeros(2,N);                                               % soft-decision phase estimate
            for k=L:N
                % Make soft decision on symbol y(k)
                yk = y(:,k);
                psik_tilde = (1/M)*angle(-yk.^M);
                p = floor(0.5+(mean(psi(:,k-3:k-1),2)-psik_tilde)/(2*pi/M));
                psi(:,k) = psik_tilde + p*(2*pi/M);
                thetahat(k-Delta) = Whd' * mean(psi(:,k:-1:k-L+1)).';
              
                % Detect symbol
                x1hat(k-Delta) = y1(k-Delta)*exp(-1i*thetahat(k-Delta));
                x2hat(k-Delta) = y2(k-Delta)*exp(-1i*thetahat(k-Delta));

            end;
        elseif ndtype == '2';
            Wsd = wiener_fir(L,Delta,sigmandsq,sigmapsq);
            % Raise received signal to M-th power and filter
            yM = y.^M;
            z = zeros(1,N);
            for k=L:N
                z(k-Delta) = Wsd' * mean(yM(:,k:-1:k-L+1)).';

                % Compute angle
                thetatilde = (1/M)*angle(-z(k-Delta));
                p = floor(0.5+(thetahat(k-Delta-1)-thetatilde)/(2*pi/M));
                thetahat(k-Delta) = thetatilde + p*(2*pi/M);
               
                % Detect symbol
                x1hat(k-Delta) = y1(k-Delta)*exp(-1i*thetahat(k-Delta));
                x2hat(k-Delta) = y2(k-Delta)*exp(-1i*thetahat(k-Delta));
            end;
        end
    end

    xhat = [x1hat.';x2hat.'];    
    
    function x = detect(y)
        x = qammod(qamdemod(y, M, 0, 'Gray'), M, 0, 'Gray');
    end
end