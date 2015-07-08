%% Transmitted power necessary to achieve target BER = 1e-4
% SOA gain was optimize for every point
% No modulator, no fiber, and receivere filter matched to pulse shape
% Opt. filter BW/Ele filter BW = 4
% Bit rate = 100 GHz
% No RIN or shot noise
% thermal noise NEP = 20 pa/sqrt(Hz)
clear, clc, close all

addpath ../../f/

M = [2, 4];
FndB = 10:-1:3;

Ptx.uniform =[[-25.3512732864707,-26.3116344172544,-27.3861483911839,-28.3024955409606,-29.2954062998511,-30.1763469411764,-31.0720478916098,-31.8975170318889];
    [-16.7411209939212,-17.6447211231479,-18.5838073052381,-19.5442311603577,-20.5201860646692,-21.5037804418418,-22.4933574962827,-23.4765192989012]];
             
Ptx.nonuniform = [[-26.9349412724686,-27.8571232188149,-28.8797169899349,-29.8069337973133,-30.7591378197701,-31.6171641040208,-32.4554965876010,-33.2691133822184];
    [-20.8703971468276,-21.8565062094265,-22.8385511294927,-23.8231503229074,-24.7995318298993,-25.7514501452422,-26.6573770354024,-27.4881357065587]];

Ptx.nonuniform_gauss_approx = [[-24.8252400249940,-26.0910530549283,-26.9986021737613,-28.0463376902151,-29.1769538180169,-30.2940609939796,-31.3531144235331,-32.3978652995545];
    [-20.2375660799052,-21.1850867366061,-22.1346630988215,-23.0934033861980,-24.0745475890736,-25.0610906742350,-26.0333398500522,-26.9743292755961]];
    
% Plot  
Colors = {'k', 'b', 'r', 'g'};

figure, hold on, box on, grid on
for k = 1:size(Ptx.uniform, 1)
    plot(FndB, Ptx.uniform(k, :), ':', 'Color', Colors{k})
    plot(FndB, Ptx.nonuniform(k, :), '-', 'Color', Colors{k})
    plot(FndB, Ptx.nonuniform_gauss_approx(k, :), '--', 'Color', Colors{k})
end

legend('Equally-spaced levels', 'Optimized level spacing', 'Optimized level spacing w/ Gaussian Approx.', 'Location', 'NorthWest')
xlabel('Noise Figure (dB)')
ylabel('Transmitted Power (dBm)')

% convert gca to latex
% matlab2tikz files must be in current directory
matlab2tikz('Ptx_x_Fn.tikz')
