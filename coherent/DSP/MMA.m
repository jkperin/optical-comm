% CMA.m    Milad Sharif  2-7-13

function [xhat, eps, WT] = MMA(tsamp, Ysq, Nsymb, oversamp, AdEq)

    Lfilt    = AdEq.LFilt;
    mu       = AdEq.mu;


    % Initialize digital filters                                                % Based on Section III P. Winzer, et all, “ Spectrally efficient 
    %WT = 0.01*ones(2,4*Lfilt+2);                                                    % long-haul optical networking using 112-Gb/s polarizationmultiplexed                       
    WT = zeros(2,4*Lfilt+2);
    WT(1,[Lfilt,Lfilt+1])=1/sqrt(2);                                            % 16-QAM,” J. Lightw. Technol., vol. 28, no. 4, pp. 547–556, Feb. 15, 2010
    WT(2,(2*Lfilt+1)+[Lfilt,Lfilt+1])=1/sqrt(2);

    % initialize filter output vectors
    x1hat = zeros(1,Nsymb);
    x2hat = zeros(1,Nsymb);


    % initialize error vectors
    % based on transmitted symbols
    eps1 = zeros(1,Nsymb);
    eps2 = zeros(1,Nsymb);

    Y1sq = Ysq(1,:);
    Y2sq = Ysq(2,:);

    % Run LMS adaptive algorithm. The index k runs over symbols.
    for k = 1:Nsymb;
        % form the Yk
        eqindk = mod(floor(oversamp*(k-1))+Lfilt:-1:floor(oversamp*(k-1))-Lfilt,size(tsamp,2))+1; % indices for Yk
        Yk = [Y1sq(eqindk),Y2sq(eqindk)].';        % form the Yk
        % select the correct WT

        % compute and store equalizer output at time k
        xkhat = WT*Yk;          % equalizer output at time k

        x1hat(k) = xkhat(1,:);  % store equalizer output 1 at time k
        x2hat(k) = xkhat(2,:);  % store equalizer output 2 at time k


        % Compute and store error at time k. Sign of error is such that when there is no ISI, error = noise
        % based on transmitted symbols

        eps1_R(k) = error(real(x1hat(k)));
        eps1_I(k) = error(imag(x1hat(k)));
        
        eps2_R(k) = error(real(x2hat(k)));
        eps2_I(k) = error(imag(x2hat(k)));
        
        % based on decision symbols
        %eps1D(k) = xkhat(1,:)-xkD(1,:);
        %eps2D(k) = xkhat(2,:)-xkD(2,:);

        % update WT
        WT = WT + mu*([eps1_R(k)*real(x1hat(k))+1j*eps1_I(k)*imag(x1hat(k))...
            ;eps2_R(k)*real(x2hat(k))+1j*eps2_R(k)*real(x2hat(k))]*Yk');
        % store WT to the correct place
    end
    % Compile outputs 
    xhat = [x1hat; x2hat];
    eps = [eps1_R+1j*eps2_I;eps2_R+1j*eps2_I];
end

function e = error(x)
    R2 = 1;
    e = R2 - abs(x)^2;
end