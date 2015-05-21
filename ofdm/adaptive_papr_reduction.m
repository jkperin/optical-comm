clear, clc, close all


Nsymb = 2^14;
Nc = 64;
Nu = 52;
CS = 64;

%% Generate OFDM signal at chip rate (done in DSP)
CS = CS*ones(1, Nc);
dataTX = zeros(Nu/2, Nsymb);
dataTXm = zeros(Nu/2, Nsymb);
for kk = 1:Nu/2
    dataTX(kk,:) = randi([0 CS(kk)-1], [1 Nsymb]);              % data to be encoded (symbols columnwise)
    dataTXm(kk,:) = qammod(dataTX(kk,:), CS(kk), 0, 'gray');    % encoded QAM symbols to be modulated onto subcarriers
%     dataTXm(kk,:) = sqrt(Pn(kk))*dataTXm(kk,:)/sqrt(Pqam(kk));  % scale constellation so that Pn is the power at nth subcarrier (E(|Xn|^2) = Pn)
end

% zero-pad and ensure Hermitian symmetry
% -f(N) ... -f(1) f(0) f(1) ... f(N-1) => (0 ... 0, X(-n)*, 0, X(n), 0 ... 0)
Xn = ifftshift([zeros((Nc-Nu)/2, Nsymb); ...     % zero-pad neg. frequencies
              flipud(conj(dataTXm));...            % data*(-n) (Nu) 
              zeros(1, Nsymb); ...                % 0 at f == 0 (1)
              dataTXm; ...                         % data(n)  (Nu)
              zeros((Nc-Nu)/2 - 1, Nsymb)], 1);   % zero-pad pos. frequencies

% % Perform ifft (Nc-IDFT) columnwise
xn = Nc*ifft(Xn, Nc, 1); 
xn = reshape(xn, 1, numel(xn));

sigtx = std(xn);

mu =0.00051;
rclip = 1;
K = 1 - 2*qfunc(rclip);
Cnp = ones(Nc/2+1, 1);
eps = zeros(Nsymb, 1);
vardn = zeros(Nsymb, 1);

for k = 2:Nsymb
    Cn = [Cnp; flipud(conj(Cnp(2:end-1)))];
    xnk = Nc*ifft(Xn(:, k).*Cn, Nc, 1);
    
    % Note: both positive and negative tails are clipped
    xnc = xnk;
    xnc(xnk < -rclip*sigtx) = -rclip*sigtx;
    xnc(xnk > rclip*sigtx) = rclip*sigtx;
    
    dn = xnc - xnk;
    
    Dn = fft(dn)/Nc;
    
    Dn = Dn(1:Nc/2+1);
    
    vardn(k) = var(dn);
    
%     eps(k) = vardn(k) - vardn(k-1);
%     eps(k) = mean(dn);
%     [m, ix] = max(abs(dn));
%     eps(k) = dn(ix);
        
    Cnp = Cnp - 2*mu*eps(k)*Xn(1:Nc/2+1, k);
    
end

plot(vardn)

figure, hold on
plot(xn, 'b')
plot([1 numel(Xn)], sigtx*rclip*[1 1], 'k')
plot([1 numel(Xn)], -sigtx*rclip*[1 1], 'k')

