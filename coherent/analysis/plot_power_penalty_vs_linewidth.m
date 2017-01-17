%% Plot power penalty vs linewidth
clear, clc, close all

addpath ../f/
addpath ../../f/

Qpsk = QAM(4, 2*112e9);
sim.ModFormat = Qpsk;
BERtarget = 1.8e-4;
sim.BERtarget = BERtarget;
totalLinewidth = 2*1e3*(0:100:1e3);
Delay = 415e-12;
csi = sqrt(2)/2;
R = 1;
q = 1.60217662e-19;

SNRdBref = SNRreq(BERtarget, Qpsk.M, Qpsk.type);
SNRdB2PrxdBm = @(SNRdB) 10*log10(10^(SNRdB/10)*2*q*Qpsk.Rs/(R*1e-3));
PrefdBm = SNRdB2PrxdBm(SNRdBref);

for Ncpr = 1:2
    for k = 1:length(totalLinewidth)
        [wnOpt(Ncpr, k), ~, SNRdBstart] = optimizePLL(csi, Delay, totalLinewidth(k), Ncpr, sim, true);
        
        [SNRdBpen, ~, exitflag] = fzero(@(SNRdB)...
            log10(ber_qpsk_imperfect_cpr(SNRdB, phase_error_variance(csi, wnOpt(Ncpr, k), Ncpr, Delay, [totalLinewidth(k), 0], SNRdB, Qpsk.Rs))) - log10(BERtarget), SNRdBstart);
        
        if exitflag ~= 1
            warning('SNRdB calculation finished with exitflag = %d', exitflag);
        end
        
        PrxreqdBm(Ncpr, k) = SNRdB2PrxdBm(SNRdBpen);
    end
end

%% BERtarget = 1.8e-4
PrxreqdBm_sim_1pol_200ps = [-33.7568543356963,-33.4144863343955,-33.3560796656632,-33.2240882178380,-32.8997299311060,-33.0926423180134,-32.9281003655659,-32.7944228772677,-32.6658656773080,-32.5536180425341,-32.5091508056068;-33.7568543356963,-33.4964021332317,-33.3531576373747,-33.2821085619065,-33.1532834975521,-33.1722278401910,-33.1089745475579,-32.9326580665056,-32.7829060179932,-32.7218532314027,-32.6537638433014];
PrxreqdBm_sim_2pol_200ps = [-33.7568543356963,-33.2787326076154,-33.2539068510371,-33.1266261226035,-33.0348971080981,-32.8024583561056,-32.5616155878530,-32.4068039018090,-32.3703265102843,-32.1467345784369,-32.0867664369898;-33.7568543356963,-33.4390632328829,-33.1268313430039,-32.9785208820970,-33.1048960758787,-32.9577380269151,-32.8443976348539,-32.7013505264218,-32.8099286269006,-32.5899658387921,-31.6759995042541];

PrxreqdBm_sim_1pol_400ps = [-33.8145225948464,-33.3656633486190,-33.2143574298784,-33.0539564044703,-32.8920365059379,-32.7927916660063,-32.4279835446653,-32.4122247434566,-31.9553078688746,-31.7920314073379,-31.4294955160951;-33.8145225948464,-33.5557413258355,-33.2844956152750,-33.1610496159534,-33.0169044455728,-32.9652482083036,-32.4617999491913,-32.5259100729796,-32.4669366027757,-32.2978904821708,-32.0625996651640]
PrxreqdBm_sim_2pol_400ps = [-33.8145225948464,-33.3807861553356,-33.2787568609794,-32.7980613758601,-32.3784723483050,-32.3563955910704,-31.8625728157591,-32.0912130592712,-31.1581751659570,-31.5131599578334,-30.8822060883747;-33.8145225948464,-33.4326273768968,-33.1491307630770,-32.9543478085719,-32.9217843542351,-32.5946555345277,-32.4129656178407,-32.0391775139557,-31.9590501384141,-31.6539265566709,-31.8338541521852];

figure(100), hold on, box on
plot(totalLinewidth/1e3, PrxreqdBm(1, :)-PrefdBm, '--k', 'LineWidth', 2)
% plot(totalLinewidth/1e3, PrxreqdBm_sim_1pol_200ps(1, :)-PrefdBm, '--ob', 'LineWidth', 2)
% plot(totalLinewidth/1e3, PrxreqdBm_sim_1pol_200ps(2, :)-PrefdBm, '--or', 'LineWidth', 2)
plot(totalLinewidth/1e3, PrxreqdBm(2, :)-PrefdBm, '-k', 'LineWidth', 2)
% plot(totalLinewidth/1e3, PrxreqdBm_sim_2pol_200ps(1, :)-PrefdBm, '-ob', 'LineWidth', 2)
% plot(totalLinewidth/1e3, PrxreqdBm_sim_2pol_200ps(2, :)-PrefdBm, '-or', 'LineWidth', 2)
% axis([totalLinewidth([1 end])/1e3 -36 -32])
xlabel('Linewidth (kHz)', 'FontSize', 12)
ylabel('Power penalty (dB)', 'FontSize', 12)
set(gca, 'FontSize', 12)
legend('Analysis', 'Simulation: XOR-based', 'Simulation: Costas') 
grid on
m2tikz = matlab2tikz(gca);
m2tikz.write('comparison_1_vs_2_pols.tikz');


% figure(101), hold on, box on
% plot(totalLinewidth/1e3, PrxreqdBm(1, :), 'k', 'LineWidth', 2)
% plot(totalLinewidth/1e3, PrxreqdBm_sim_1pol_400ps(1, :), '-b', 'LineWidth', 2)
% plot(totalLinewidth/1e3, PrxreqdBm_sim_1pol_400ps(2, :), '-r', 'LineWidth', 2)
% plot(totalLinewidth/1e3, PrxreqdBm(2, :), '--k', 'LineWidth', 2)
% plot(totalLinewidth/1e3, PrxreqdBm_sim_2pol_400ps(1, :), '--b', 'LineWidth', 2)
% plot(totalLinewidth/1e3, PrxreqdBm_sim_2pol_400ps(2, :), '--r', 'LineWidth', 2)
% axis([totalLinewidth([1 end])/1e3 -36 -32])
% xlabel('Linewidth (kHz)', 'FontSize', 12)
% ylabel('Receiver sensitivity (dBm)', 'FontSize', 12)
% set(gca, 'FontSize', 12)
% legend('Analysis', 'Simulation: XOR-based', 'Simulation: Costas') 
% grid on
% 
% figure, hold on, box on
% plot(totalLinewidth, wnOpt/(2*pi*1e6))
