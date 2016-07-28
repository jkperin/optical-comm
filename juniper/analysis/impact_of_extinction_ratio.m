%% Impact of extinction ratio

addpath ../../mpam/
addpath ../f/

mpam = PAM(4, 56e9);

OSNRdB = linspace(25, 35);

ERdB = [-5 -8 -10 -12 -15 -20 -30 -40 -Inf];

figure, hold on, box on

leg = {};
for k = 1:length(ERdB)
    plot(OSNRdB, log10(pam_ber_from_osnr(mpam.M, OSNRdB, mpam.Rs/2, -ERdB(k))), 'linewidth', 2)
    leg = [leg sprintf('ER = %d dB', -ERdB(k))];
end
OSNRdBexp = [25,26,27,28,29,30,31,32,34,36,39,43,48];
BERcountexp = [0.00236459514949504,0.00134356948798948,0.000753756517317816,0.000367979274457639,0.000183607356782545,0.000138263029803860,4.52555156087015e-05,2.67081731461189e-05,5.19325588952313e-06,2.97339849040559e-06,7.41893698503304e-07,0,0];
plot(OSNRdBexp, log10(BERcountexp), '-ok', 'linewidth', 2)
leg = [leg 'Experiment'];
xlabel('OSNR (dB)')
ylabel('log_{10} (BER)')
legend(leg)
axis([25 35 -8 0])
