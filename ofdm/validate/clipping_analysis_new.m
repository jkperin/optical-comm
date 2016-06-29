%% Validate expressions for the 1st and 2nd moments of a clipped normal-distributed signal
clear, clc, close all

N = 2^16;       % number of points
s = 2.5476;     % std

%% Clipping once
r = 0.5;           % clipping ratio

for k = 1:100
    x = s*randn(1, N);
    xc = x;
    xc(x < -r*s) = -r*s;
    d = xc - x*(1-qfunc(r));
    xc = xc + r*s;
    Exc(k) = mean(xc);
    Exc2(k) = mean(abs(xc).^2);
    varxc(k) = var(xc);
    
    Ed(k) = mean(d);
    Ed2(k) = mean(abs(d).^2);
    vard(k) = var(d);
end

disp('Clipping once')
[mean(Exc), s*(r*(1-qfunc(r)) + 1/sqrt(2*pi)*exp(-r^2/2))]
[mean(Exc2), s^2*((1+r^2)*(1-qfunc(r))+ r/sqrt(2*pi)*exp(-r^2/2))]
[mean(varxc), s^2*((1-qfunc(r))*(1 + r^2*qfunc(r)) - r/sqrt(2*pi)*exp(-r^2/2)*(1-2*qfunc(r)) - 1/(2*pi)*exp(-r^2))]


[mean(Ed), s*(-r*qfunc(r) + 1/sqrt(2*pi)*exp(-r^2/2))]
[mean(Ed2), s^2*(qfunc(r)*(1 + r^2) - qfunc(r)^2 - r/sqrt(2*pi)*exp(-r^2/2))]
[mean(vard), s^2*(qfunc(r)*(1-qfunc(r))*(r^2+1) - r/sqrt(2*pi)*exp(-r^2/2)*(1-2*qfunc(r)) - 1/(2*pi)*exp(-r^2))]

% [mean(Exc), s*(-r*qfunc(r) + 1/sqrt(2*pi)*exp(-r^2/2))]
% [mean(Exc2), s^2*(r^2*qfunc(r) + 1 - qfunc(r) - r/sqrt(2*pi)*exp(-r^2/2))]

clear Ex* Ed* var*

%% Clipping twice (DC-OFDM)
r1 = 1;           % clipping ratio
r2 = 1;

K = 1 - qfunc(r1) - qfunc(r2);

for k = 1:100
    rng('shuffle')
    x = s*randn(1, N);
    xc = x;
    xc(x < -r1*s) = -r1*s;
    xc(x > r2*s) = r2*s;
    d = xc - x*K;
    
    % Note that dc bias was not added. Expressions are written for this
    % case
    Exc(k) = mean(xc);
    Exc2(k) = mean(abs(xc).^2);
    varxc(k) = var(xc);
    
    Ed(k) = mean(d);
    Ed2(k) = mean(abs(d).^2);
    vard(k) = var(d);
end

disp('Clipping twice')
[mean(Exc), s*(-r1*qfunc(r1) + r2*qfunc(r2) + 1/sqrt(2*pi)*(exp(-r1^2/2) - exp(-r2^2/2)))]
[mean(Exc2), s^2*(1 + qfunc(r1)*(r1^2 - 1) + qfunc(r2)*(r2^2 - 1) - 1/sqrt(2*pi)*(r1*exp(-r1^2/2) + r2*exp(-r2^2/2)))]
% only for r1 = r2
[mean(varxc), s^2*(1 + 2*qfunc(r1)*(r1^2 - 1) - 2*r1/sqrt(2*pi)*exp(-r1^2/2))]


[mean(vard), s^2*(1 - K^2 + 2*qfunc(r1)*(r1^2 - 1) - 2*r1/sqrt(2*pi)*exp(-r1^2/2))]


% [mean(Ed), s*(-r*qfunc(r) + 1/sqrt(2*pi)*exp(-r^2/2))]
% [mean(Ed2), s^2*(qfunc(r)*(1 + r^2) - qfunc(r)^2 - r/sqrt(2*pi)*exp(-r^2/2))]
% [mean(vard), s^2*(qfunc(r)*(r^2+1) + qfunc(r)^2*(r^2-1) - r/sqrt(2*pi)*exp(-r^2/2)*(1-2*qfunc(r)) - 1/(2*pi)*exp(-r^2))]

% [mean(Exc), s*(-r*qfunc(r) + 1/sqrt(2*pi)*exp(-r^2/2))]
% [mean(Exc2), s^2*(r^2*qfunc(r) + 1 - qfunc(r) - r/sqrt(2*pi)*exp(-r^2/2))]


%% Clipping twice (ACO-OFDM)
disp('Clipping twice (ACO-OFDM)');
r1 = 0;           % clipping ratio
r2 = 1;

K = 1 - qfunc(r1) - qfunc(r2);

for k = 1:100
    rng('shuffle')
    x = s*randn(1, N);
    xc = x;
    xc(x < -r1*s) = -r1*s;
    xc(x > r2*s) = r2*s;
    d = xc - x*K;
    
    % Note that dc bias was not added. Expressions are written for this
    % case
    Exc(k) = mean(xc);
    Exc2(k) = mean(abs(xc).^2);
    varxc(k) = var(xc);
    
    Ed(k) = mean(d);
    Ed2(k) = mean(abs(d).^2);
    vard(k) = var(d);
end

[mean(Exc), s*(r2*qfunc(r2) + 1/sqrt(2*pi)*(1 - exp(-r2^2/2)))]
[mean(Exc2), s^2*(1/2 + qfunc(r2)*(r2^2 - 1) - r2/sqrt(2*pi)*exp(-r2^2/2))]

[mean(vard), s^2*(1/2 + qfunc(r2)*(r2^2 - 1) - r2/sqrt(2*pi)*exp(-r2^2/2) - (r2*qfunc(r2) + 1/sqrt(2*pi)*(1 - exp(-r2^2/2)))^2 - K^2)]


% [mean(Ed), s*(-r*qfunc(r) + 1/sqrt(2*pi)*exp(-r^2/2))]
% [mean(Ed2), s^2*(qfunc(r)*(1 + r^2) - qfunc(r)^2 - r/sqrt(2*pi)*exp(-r^2/2))]
% [mean(vard), s^2*(qfunc(r)*(r^2+1) + qfunc(r)^2*(r^2-1) - r/sqrt(2*pi)*exp(-r^2/2)*(1-2*qfunc(r)) - 1/(2*pi)*exp(-r^2))]

% [mean(Exc), s*(-r*qfunc(r) + 1/sqrt(2*pi)*exp(-r^2/2))]
% [mean(Exc2), s^2*(r^2*qfunc(r) + 1 - qfunc(r) - r/sqrt(2*pi)*exp(-r^2/2))]





% r = linspace(0, 6);

% figure
% plot(r, r.*(1-qfunc(r)) + 1/sqrt(2*pi)*exp(-r.^2/2), 'k', 'LineWidth', 2)
% hold on
% plot(r, -r.*qfunc(r) + 1/sqrt(2*pi)*exp(-r.^2/2), 'r', 'LineWidth', 2)
% xlabel('Clipping ratio (r)', 'FontSize', 12)
% ylabel('Normalized Mean', 'FontSize', 12)
% legend('Clipped signal', 'Clipping noise')
% 
% 
% figure
% plot(r, (1+r.^2).*(1-qfunc(r))+ r/sqrt(2*pi).*exp(-r.^2/2), 'k', 'LineWidth', 2)
% hold on
% plot(r, qfunc(r).*(1 + r.^2) - qfunc(r).^2 - r/sqrt(2*pi).*exp(-r.^2/2), 'r', 'LineWidth', 2)
% xlabel('Clipping ratio (r)', 'FontSize', 12)
% ylabel('Normalized Power', 'FontSize', 12)
% legend('Clipped signal', 'Clipping noise')

% figure
% plot(r, (1-qfunc(r)).*(1 + r.^2.*qfunc(r)) - r/sqrt(2*pi).*exp(-r.^2/2).*(1-2*qfunc(r)) - 1/(2*pi)*exp(-r.^2), 'k', 'LineWidth', 2)
% hold on
% plot(r, qfunc(r).*(r.^2+1) + qfunc(r).^2.*(r.^2-1) - r/sqrt(2*pi).*exp(-r.^2/2).*(1-2*qfunc(r)) - 1/(2*pi)*exp(-r.^2), 'r', 'LineWidth', 2)
% xlabel('Clipping ratio (r)', 'FontSize', 12)
% ylabel('Normalized Variance (dB)', 'FontSize', 12)
% legend('Clipped signal', 'Clipping noise')
% 
% figure
% plot(r, 10*log10((1-qfunc(r)).^2./(qfunc(r).*(r.^2+1) + qfunc(r).^2.*(r.^2-1) - r/sqrt(2*pi).*exp(-r.^2/2).*(1-2*qfunc(r)) - 1/(2*pi)*exp(-r.^2))), 'k', 'LineWidth', 2)
% xlabel('Clipping ratio (r)', 'FontSize', 12)
% ylabel('SNR total (dB)', 'FontSize', 12)
% legend('Clipped signal', 'Clipping noise')
% 
% NEP = 30e-12;
% N0 = NEP^2/2;
% fs = (40:20:120)*1e9;
% ros = 1.23;
% Nc = 64;
% Nu = 52;
% sig2 = (2.0972e-04)^2;
% varth = ((N0*fs/ros)/Nc)/sig2;
% 
% figure
% plot(r, 10*log10((1/Nu*(1-qfunc(r)).^2)./...
% (varth(end)+ 1/Nc*(qfunc(r).*(r.^2+1) + qfunc(r).^2.*(r.^2-1) - r/sqrt(2*pi).*exp(-r.^2/2).*(1-2*qfunc(r)) - 1/(2*pi)*exp(-r.^2)))), 'k', 'LineWidth', 2)
% hold on
% plot(r, 17.8504*ones(size(r)), '--b') % SNR to get to 1.8e-4 using 16-QAM
% plot(r, 23.9128*ones(size(r)), '--c') % SNR to get to 1.8e-4 using 16-QAM
% xlabel('Clipping ratio (r)', 'FontSize', 12)
% ylabel('SNR (dB)', 'FontSize', 12)
