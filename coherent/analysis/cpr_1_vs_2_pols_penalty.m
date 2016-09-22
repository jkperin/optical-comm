%
clear, clc

addpath ../f/
addpath ../../f/

Qpsk = QAM(4, 2*112e9);
sim.ModFormat = Qpsk;
BERtarget = 1e-3;
sim.BERtarget = BERtarget;
totalLinewidth = 2e3*(0:200:2e3);
Delay = 200e-12;
csi = sqrt(2)/2;
q = 1.60217662e-19;

for Ncpr = 1:2
    for k = 1:length(totalLinewidth)
        [wnOpt(Ncpr, k), ~, phiError, SNRdB] = optimizePLL(csi, Delay, totalLinewidth(k), Ncpr, sim, true);
        
        % PrxrefdBm = 10*log10(10^(SNRdB/10)*2*q*Qpsk.Rs/1e-3)

        [SNRdBpen, ~, exitflag] = fzero(@(SNRdB)...
            log10(ber_qpsk_imperfect_cpr(SNRdB, phase_error_variance(csi, wnOpt(Ncpr, k), Ncpr, Delay, totalLinewidth(k), SNRdB, Qpsk.Rs))) - log10(BERtarget), SNRdB);
        
        if exitflag ~= 1
            warning('SNRdB calculation finished with exitflag = %d', exitflag);
        end
        
        PrxreqdBm(Ncpr, k) = 10*log10(10^(SNRdBpen/10)*2*q*Qpsk.Rs/1e-3);
    end
end

PrxreqdBm_sim = [-35.1345468969015,-34.7499210561243,-34.5445802999749,-34.4030302511419,-34.3236275896456,-34.4503660505951,-23.3331085328622,-34.0525455391051,-33.6408219836450,-33.4006842227610,-33.1347171804686;-35.1345468969015,-34.8174775394835,-34.6844633477978,-34.5912964234285,-34.5766182352072,-33.9533188162426,-34.3866039942583,-34.2906561271326,-34.0419234003847,-33.9211750295863,-33.7536301259386];

figure(100), hold on, box on
plot(totalLinewidth, PrxreqdBm)
plot(totalLinewidth, PrxreqdBm_sim)
axis([totalLinewidth([1 end]) -36 -32])
        
figure, hold on, box on
plot(totalLinewidth, wnOpt/(2*pi*1e6))
