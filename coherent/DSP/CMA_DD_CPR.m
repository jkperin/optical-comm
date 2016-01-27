% adeq.m    Milad Sharif  6-25-13
% Time-domain linear equalizer adapted by LMS algorithm.
% Only for  mod(oversamp*k,1) == 0 for every k
function [xhat, xD, eps, WT, thetahat, Delta] = CMA_DD_CPR(tsamp, Ysq, Nsymb, oversamp, AdEq,Cpr , PN, snr_dB)

    Lfilt    = AdEq.LFilt;

    muVec    = AdEq.muVec;
    hardSwitch = AdEq.hardSwitch;
    
    % initialize digital filters ------------------------------------------
    WT = zeros(2,4*Lfilt+2);
    WT(1,[Lfilt,Lfilt+1])=.1/sqrt(2);                                            % 16-QAM,” J. Lightw. Technol., vol. 28, no. 4, pp. 547–556, Feb. 15, 2010
    WT(2,(2*Lfilt+1)+[Lfilt,Lfilt+1])=.1/sqrt(2);

    % initialize filter output vectors ------------------------------------
    xequ = zeros(2,Nsymb);
    xhat = zeros(2,Nsymb);
    
    % initialize decision vectors -----------------------------------------
    xD = zeros(2,Nsymb);

    % initialize error vectors based on decision symbols ------------------
    eps = zeros(2,Nsymb);

    % split into I and Q
    Y1sq = Ysq(1,:);
    Y2sq = Ysq(2,:);

    % Carrier phase recovery----------------------------------------------- 
    psi = zeros(2,Nsymb);                                                       % soft-decision phase estimate
    theta_tilde = zeros(2,Nsymb);                                               % input to soft-decision phase estimator
    thetahat = zeros(Nsymb,1);                                                  % hard-decision phase estimate
        
    sigmapsq = PN.sigma_p^2;
    crtype   = Cpr.type;
    Lmax     = Cpr.Filter.Lmax;
    fract    = Cpr.Filter.f;
    filttype = Cpr.Filter.Type;

    M = 4;                                                                  % M-PSK
    Ptx = 1;                                                                % transmitted power in one polarization

    % perform some computations related to carrier recovery
    a = sqrt(Ptx/2);                                                        % symbol spacing = 2a
    absx = abs(a*[1+1i*1]);
    eta_c = sum(absx.^2) * sum(absx.^(-2));                                 % constellation penalty

    snr_des = 10^(0.1*snr_dB);                                              % SNR per symbol at which the carrier recovery is optimized

    if strcmp(crtype,'DD')
        sigmandsq = 0.5*eta_c/snr_des;                                      % variance of n'_k
    elseif strcmp(crtype,'NDA')
        sigmandsq = (1/(2*M^2))*sum(((gamma(M+1)./(gamma([1:M]+1).*gamma(M-[1:M]+1))).^2).*gamma([1:M]+1)./(snr_des.^[1:M]));
    end

    r = sigmapsq/sigmandsq;
    alpha = (1+0.5*r)-sqrt((1+0.5*r)^2-1);                                  % optimum IIR coefficient

    if strcmp(filttype,'FIR')
        Lfract = ceil(-2*log(1-fract)/(log(1/alpha)));                      % required length of FIR filter to include a fraction fract of the total energy
        L = min(Lmax,Lfract);                                               % actual length of Whd (for DD or NDA 1) or of Wsd (for NDA 2) employed
        Delta = floor(0.5*(L-1));                                           % delay of Whd (for DD or NDA 1) or of Wsd (for NDA 2)
        if strcmp(crtype,'DD')
            d = floor(0.5*(L-1));                                           % length of Wsd for DD
        end
    elseif strcmp(filttype,'IIR')
        Lfract = ceil(-2*log(1-fract)/(log(1/alpha)));                      % required length of FIR filter to include a fraction fract of the total energy
        L = Lfract;                                                         % length of Whd (for DD or NDA 1) or of Wsd (for NDA 2) employed in approximating the IIR filter
        Delta = 0;                           
        if strcmp(crtype,'DD')
            d = floor(0.5*((min(Lfract,Lmax)-1)));                          % length of Wsd for DD
        end
    end
    
    Wsd = wiener_fir(d,0,sigmandsq,sigmapsq);
    Whd = wiener_fir(L,Delta,sigmandsq,sigmapsq);       
    
    % Run LMS adaptive algorithm. The index k runs over symbols. ----------
    for k = 1:Nsymb;
        
%         imagesc(abs(WT));
%         colormap hot; axis image
%         title(['k: ' num2str(k) ', equalizer coefficients |\itW\rm^{T}|'])
%         pause(.05)
        
        % form the Yk 
        eqindk = mod(floor(oversamp*(k-1))+Lfilt:-1:floor(oversamp*(k-1))-Lfilt,size(tsamp,2))+1;
        Yk = [Y1sq(eqindk),Y2sq(eqindk)].';                                 % form the Yk
        % compute and store equalizer output at time k
        xkequ = WT*Yk;                                                      % equalizer output at time k    
        xequ(:,k) = xkequ;                                                  % store equalizer output at time k

        % Make soft decision on symbol xequ(k) ============================
        if (k>=L)
            sk = xkequ.*exp(-1i*theta_tilde(:,k));
            xhatk_sd   = detect(sk);
            psik_tilde = angle(xkequ)-angle(xhatk_sd);
            p  = floor((psi(:,k-1)-psik_tilde+pi)/(2*pi));
            psi(:,k) = psik_tilde + p*(2*pi);
            theta_tilde(:,k+1) = psi(:,k:-1:k-d+1) * Wsd;
            thetahat(k-Delta) = Whd' * mean(psi(:,k:-1:k-L+1)).';

            ind = k-Delta;
            xhat(:,ind) = xequ(:,ind)*exp(-1i*thetahat(ind));
            xD(:,ind) = detect(xhat(:,ind));                        % decision at time k            
            
            % based on decision symbols
            eps(:,ind) = xhat(:,ind) - xD(:,ind);
        end
        
        
        % update WT =======================================================
        muk = muVec(k);
        if (k < hardSwitch)
            epsk = 2 - abs(xkequ).^2;
            eps(:,k) = epsk;
            WT = WT + muk*([epsk(1)*xkequ(1)*Yk';epsk(2)*xkequ(2)*Yk']);
            
        elseif mod(k-L,Delta+1) == 0
            ind = k-Delta;
            eqind = mod(floor(oversamp*(ind-1))+Lfilt:-1:floor(oversamp*(ind-1))-Lfilt,size(tsamp,2))+1;
            Yind = [Y1sq(eqind),Y2sq(eqind)].';
            % step size parameter at time k
            % error at time k (based on transmitted symbols or decision symbols, depending on value of decwtvec
            epsk = eps(:,ind);            
            % update WT
            WT = WT - 2*muk*([epsk(1)*Yind';epsk(2)*Yind']);
            % store WT to the correct place
        end
    end
end
function x = detect(y)
    x = sign(real(y))+1i*sign(imag(y));
end

