%% test apd class

clear, clc, close all

GainBW = 390e9;
ka = 0.09;
G = 10;
dt = 1/50e9;

rxg = apd(GainBW, ka, G, 'Gaussian');
rxd = apd(GainBW, ka, G, 'DoublyStoch');

PxdBm = -30:2:-20;
Px = 1e-3*10.^(PxdBm/10);
figure, hold on
for k = 1:length(Px)
    lambda = (rxd.R*Px(k) + rxd.Id)*(dt/2)/rxd.q;

    [pk_sp, n1] = rxd.output_pmf_saddlepoint(lambda);
%     pk_ex = rxd.output_pmf(lambda);
    n1 = n1*rxd.q/(dt/2);

    var_shot = 2*rxd.q*rxd.Gain^2*rxd.Fa*(Px(k))*1/dt;
    px_g = pdf('normal', n1, Px(k)*rxd.Gain, sqrt(var_shot));
    
    figure(1), hold on
    plot(n1*1e3, (dt/2)/rxd.q*pk_sp, '-k', n1*1e3, px_g, '--r');
    figure(2), hold on
    plot(n1*1e3, log((dt/2)/rxd.q*pk_sp), '-k', n1*1e3, log(px_g), '--r');
    figure(3), hold on
    plot(n1*1e3, cumtrapz(n1, (dt/2)/rxd.q*pk_sp), '-k', n1*1e3, cumtrapz(n1, px_g), '--r');
    1;
end

figure(1)
legend('Saddlepoint approximation', 'Gaussian')
xlabel('Current (mA)', 'FontSize', 12)
ylabel('pdf', 'FontSize', 12)
figure(2)
legend('Saddlepoint approximation', 'Gaussian', 'Location', 'SouthWest')
xlabel('Current (mA)', 'FontSize', 12)
ylabel('log(pdf)', 'FontSize', 12)
axis([0 0.15 -20 15])

% pk = rxd.gain_distribution(lambda);
% 
% figure, hold on
% plot(0:length(pk)-1, pk, 'k')
% plot(n, pk_sp)
% 
% n = 0:length(pk_sp)-1;
% mean = sum(pk_sp.*n)
% var = sum(pk_sp.*(n-lambda*rx.Gain).^2)


% Pin = 1e-4*(pammod(randi([0, 3], 2^14, 1), 4) + 3);
% Nsymb = 
Pin = 1e-6*ones(1, Nsymb);
dt = 1/50e9;

output_g = rxg.detect(Pin, dt);
output_d = rxd.detect(Pin, dt);

figure, hist(output_g, 50)
figure, hist(output_d, 50)
