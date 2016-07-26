%% Impact of MZM nonlinearity

addpath ../mpam/
addpath f/

clear, clc

OSNRdBexp = [25,26,27,28,29,30,31,32,34,36,39,43,48];
BERcountexp = [0.00236459514949504,0.00134356948798948,0.000753756517317816,0.000367979274457639,0.000183607356782545,0.000138263029803860,4.52555156087015e-05,2.67081731461189e-05,5.19325588952313e-06,2.97339849040559e-06,7.41893698503304e-07,0,0];

OSNRdB = linspace(OSNRdBexp(1), OSNRdBexp(end));
BERtheory = pam_ber_from_osnr_mzmnl(4, OSNRdB, 56e9/4, [0.31 0.5]);
BERideal = pam_ber_from_osnr(4, OSNRdB, 56e9/4);

figure(1), hold on, box on
plot(OSNRdB, log10(BERideal), 'linewidth', 2)
plot(OSNRdB, log10(BERtheory), 'linewidth', 2)
plot(OSNRdBexp, log10(BERcountexp), '-ok', 'linewidth', 2)
% leg = [leg 'Experiment'];
% legend(leg)
xlabel('OSNR (dB)')
ylabel('log_{10} (BER)')
axis([25 40 -8 0])