function [Y, W, MSE] = lms_td_fse(X, eq, verbose)
%% Adaptive time-domain fractionally spaced linear equalizer adapted using LMS algorithm
% Inputs:
% - eq : equalizer parameters {eq.ros = oversampling ratio, eq.Ntrain = 
% number of traning symbols, eq.mu = adaptation rate}
% - X : input samples in the two polarizations [2 x N] at rate ros x Rs
% - Xtrain : symbols for training length(Xtrain) >= Ntrain
% - verbose (optional, default=false): whether to plot equalizer frequency
% response and convergence
% Output:
% - Y : Equalized symbols at rate 1 x Rs
% - Wx, Wy : Filters coefficients 
% - MSE: mean square error
M = eq.M;  % qam order
Ytrain = eq.trainSeq;
ros = eq.ros;
Ntrain = eq.Ntrain;
mu = eq.mu;

if mod(eq.Ntaps, 2) == 0 % always uses odd number of taps
    eq.Ntaps = eq.Ntaps + 1;
end
Ntaps = eq.Ntaps;

% Number of filters
[Nsamp, Nfilters]  = rat(ros);

if strcmpi(eq.structure, '4 filters')
    % Initialize filters coefficients
    Wx = zeros(2*Ntaps, Nfilters); % used for x polarization
    Wy = zeros(2*Ntaps, Nfilters); % used for y polarization

    % Initialize filters as if there were no ISI and there was no mixing between the two pols
    Wx(ceil(Ntaps/2), :) = 1; 
    Wy(ceil(Ntaps/2)+Ntaps, :) = 1;

    Nsymb = length(Ytrain);
    ex = zeros(1, Nsymb);
    ey = zeros(1, Nsymb);
    yx_hat = zeros(1, Nsymb); % output in x pol stream
    yy_hat = zeros(1, Nsymb); % output in y pol stream
    filt = 1;
    window = (1:Ntaps) - ceil(Ntaps/2); % {-(Ntaps-1)/2, ..., 0, ..., (Ntaps-1)/2}.
    kstart = ceil(Ntaps/(2*Nsamp))*Nsamp+1; 
    k = kstart; % runs over samples
    for n = ((kstart-1)/ros+1):Nsymb-ceil(Ntaps/2) % n runs over symbols
        % k indexes samples and n indexes symbols
        z = [X(1, k + window), X(2, k + window)]; % input of filter [x pol, y pol]

        % Filter output
        yx_hat(n) = z*Wx(:, filt);
        yy_hat(n) = z*Wy(:, filt);

        % Calculate error
        if n < Ntrain % Training based on training sequence
            ex(n) = yx_hat(n) - Ytrain(1, n);
            ey(n) = yy_hat(n) - Ytrain(2, n);
        else % based on deteced symbols
            ex(n) = yx_hat(n) - qammod(qamdemod(yx_hat(n), M, 0, 'Gray'), M, 0, 'Gray');
            ey(n) = yy_hat(n) - qammod(qamdemod(yy_hat(n), M, 0, 'Gray'), M, 0, 'Gray');
        end

        % Update filters coefficients
        Wx(:, filt) = Wx(:, filt) - 2*mu*ex(n)*z';
        Wy(:, filt) = Wy(:, filt) - 2*mu*ey(n)*z';

        % Increment filter index
        filt = mod(filt, Nfilters) + 1;

        % Controls how we run over the index k (whether to increment by 1 or 2)
        % This ensures that the center of the window of samples remains 
        % close to the nth symbol
        % this is periodic: e.g., ros = 5/4 -> ..., 1, 1, 2, 1,...
        if abs((k+1)/ros - n-1) > 0.5
            k = k + 2;
        else
            k = k + 1;
        end
    end
    
    W = {Wx, Wy};
    
elseif strcmpi(eq.structure, '2 filters')
    % Initialize filters coefficients
    Wx = zeros(Ntaps, Nfilters); % used for x polarization
    Wy = zeros(Ntaps, Nfilters); % used for y polarization
    Wmix = zeros(2, Nfilters);  % crossover filters

    % Initialize filters as if there were no ISI and there was no mixing between the two pols
    Wx(ceil(Ntaps/2), :) = 1; 
    Wy(ceil(Ntaps/2), :) = 1;


    Nsymb = length(Ytrain);
    ex = zeros(1, Nsymb);
    ey = zeros(1, Nsymb);
    yx_hat = zeros(1, Nsymb); % output in x pol stream
    yy_hat = zeros(1, Nsymb); % output in y pol stream
    filt = 1;
    window = (1:Ntaps) - ceil(Ntaps/2); % {-(Ntaps-1)/2, ..., 0, ..., (Ntaps-1)/2}.
    kstart = ceil(Ntaps/(2*Nsamp))*Nsamp+1; 
    k = kstart; % runs over samples
    for n = ((kstart-1)/ros+1):Nsymb-ceil(Ntaps/2) % n runs over symbols
        % k indexes samples and n indexes symbols
        zx = X(1, k + window);
        zy = X(2, k + window); 

        % Filter output
        yx = zx*Wx(:, filt);
        yy = zy*Wy(:, filt);

        yx_hat(n) = yx + yy*Wmix(1, filt);
        yy_hat(n) = yx*Wmix(2, filt) + yy;

        % Calculate error
        if n < Ntrain % Training based on training sequence
            ex(n) = yx_hat(n) - Ytrain(1, n);
            ey(n) = yy_hat(n) - Ytrain(2, n);        
        else % based on deteced symbols
            ex(n) = yx_hat(n) - qammod(qamdemod(yx_hat(n), M, 0, 'Gray'), M, 0, 'Gray');
            ey(n) = yy_hat(n) - qammod(qamdemod(yy_hat(n), M, 0, 'Gray'), M, 0, 'Gray');
        end

        % Update filters coefficients
    %     ind = ceil(Ntaps/2);
    %     tmpx = Wx(ind, filt);
    %     tmpy = Wy(ind, filt);
        Wx(:, filt) = Wx(:, filt) - 2*mu*ex(n)*zx';
        Wy(:, filt) = Wy(:, filt) - 2*mu*ey(n)*zy';
        Wmix(1, filt) = Wmix(1, filt) - 2*mu*ex(n)*conj(yy);    
        Wmix(2, filt) = Wmix(2, filt) - 2*mu*ey(n)*conj(yx);    

    %     Wmix(1, filt) = mean((Wy(ind, filt)*Wmix(1, filt) - 2*mu*ex(n)*zy(ind)')./tmpy);
    %     Wmix(2, filt) = mean((Wx(ind, filt)*Wmix(2, filt) - 2*mu*ey(n)*zx(ind)')./tmpx);
    %             
        % Increment filter index
        filt = mod(filt, Nfilters) + 1;

        % Controls how we run over the index k (whether to increment by 1 or 2)
        % This ensures that the center of the window of samples remains 
        % close to the nth symbol
        % this is periodic: e.g., ros = 5/4 -> ..., 1, 1, 2, 1,...
        if abs((k+1)/ros - n-1) > 0.5
            k = k + 2;
        else
            k = k + 1;
        end
    end
    
    W = {Wx, Wy, Wmix};
else
    error('lms_td_fse/structure not defined')
end

% Build outputs
Y = [yx_hat; yy_hat];
MSE = [abs(ex).^2; abs(ey).^2];


% Equalizer frequency response and convergence-----------------------------
if exist('verbose', 'var') && verbose
    Wx = W{1};
    Wy = W{2};
    figure(102)
    subplot(221), hold on, box on
    plot(1:length(MSE), MSE(1, :))
    a = axis;
    plot([1 1]*eq.Ntrain, [a(3) a(4)], '--k')
    axis([1 length(MSE) a(3:4)])
    xlabel('Iteration'); ylabel('MSE')
    
    subplot(222), hold on, box on
    plot(1:length(MSE), MSE(2, :))
    a = axis;
    plot([1 1]*eq.Ntrain, [a(3) a(4)], '--k')
    axis([1 length(MSE) a(3:4)])
    xlabel('Iteration'); ylabel('MSE')
    
    subplot(223), hold on, box on
    for filt = 1:size(Wx, 2)
        [hx, w] = freqz(Wx(:, filt), 1);
        hy = freqz(Wy(:, filt), 1, w);
        hline(1) = plot(w/(2*pi), abs(hx).^2, '-');
        hline(2) = plot(w/(2*pi), abs(hy).^2, '--', 'Color', get(hline(1), 'Color'));
    end
    xlabel('Frequency')
    ylabel('Magnitude')
    legend(hline, {'X pol', 'Y pol'})
    title(sprintf('%s, %d taps', eq.type, eq.Ntaps))
    
    subplot(224), hold on, box on
    for filt = 1:size(Wx, 2)
        [hx, w] = freqz(Wx(:, filt), 1);
        hy = freqz(Wy(:, filt), 1, w);
        hline(1) = plot(w/(2*pi),  unwrap(angle(hx)), '-');
        hline(2) = plot(w/(2*pi),  unwrap(angle(hy)), '--', 'Color', get(hline(1), 'Color'));
    end    
    xlabel('Frequency')
    ylabel('Phase')
    legend('X pol', 'Y pol')
    title(sprintf('%s, %d taps', eq.type, eq.Ntaps))
end
