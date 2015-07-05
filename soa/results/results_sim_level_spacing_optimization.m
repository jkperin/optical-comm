%% Transmitted power necessary to achieve target BER = 1e-4
% SOA gain was optimize for every point
% No modulator, no fiber, and receivere filter matched to pulse shape
% Opt. filter BW/Ele filter BW = 4
% Bit rate = 100 GHz
% No RIN or shot noise
% thermal noise NEP = 20 pa/sqrt(Hz)
clear, clc, close all

M = [2, 4];
FndB = 3:10;

Gsoa_2PAM_uniform = [20, 20, 20, 20, 20, 20, 20, 20];
Gsoa_2PAM_nonuniform = [20, 20, 20, 20, 20, 20, 20, 20];

Gsoa_4PAM_uniform = [20, 20, 20, 20, 20, 19.3219, 17.3356, 15.3573];
Gsoa_4PAM_nonuniform = [20, 20, 20, 20, 20, 20, 20, 20];

Ptx_nonuniform =[-35.053842044722522 -29.689049206589171;
                 -34.530321773847277 -29.172423367881063;
                 -33.951559927334543 -28.382576592343824;
                 -33.278188516684800 -27.581341770656504;
                 -32.461847978552349 -26.766494476386164;
                 -31.621100941522354 -25.870959220207215;
                 -30.756749681354776 -24.923271511066481;
                 -29.810103848061505 -23.947644063571062];
             
Ptx_uniform =   [-33.975726427211875 -26.619115350548167;
                 -33.347201126998058 -25.699381765628001;
                 -32.643459866718842 -24.740153967633447;
                 -31.888547976363803 -23.763745232141247;
                 -31.058623443278371 -22.781110390077433;
                 -30.157996521924538 -21.792366640707638;
                 -29.270686304593266 -20.807518323141053;
                 -28.280052898618191 -19.830485107051079];
             
% Plot  
Colors = {'k', 'b', 'r', 'g'};

figure, hold on, box on, grid on
for k = 1:size(Ptx_uniform, 2)
    plot(FndB, Ptx_uniform(:, k), '-', 'Color', Colors{k})
    plot(FndB, Ptx_nonuniform(:, k), '--', 'Color', Colors{k})
end

legend('Uniform Level Spacing', 'Non-Uniform Level Spacing', 'Location', 'NorthWest')
xlabel('Noise Figure (dB)')
ylabel('Transmitted Power (dBm)')

%
figure, hold on, box on, grid on
leg = cell(size(M));
for k = 1:size(Ptx_uniform, 2)
    plot(FndB, Ptx_uniform(:, k) - Ptx_nonuniform(:, k), '-', 'Color', Colors{k})
    leg{k} = sprintf('M = %d', M(k));
end

legend(leg, 'Location', 'NorthWest')
xlabel('Noise Figure (dB)')
ylabel('Power Margin Improvement (dB)')



