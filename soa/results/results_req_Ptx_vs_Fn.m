clear, clc, close all


load results_4PAM
results(1).GsoadB_opt = GsoadB_opt;
results(1).Preq = Preq;
results(1).Preq20dB.eq_spaced = [-23.676130924700775,-22.699918565227463,-21.716125208389577,-20.727182311905416,-19.734815620937407,-18.740137553107240,-17.743893199585024,-16.746577031492260];
results(1).Preq20dB.optimized = [-26.605266004025670,-25.664331557038680,-24.712016814924663,-23.744065906106474,-22.764632916790372,-21.777826671313260,-20.786886404535860,-19.792957911685612];
results(1).Preq15dB.eq_spaced = [-23.320904484745906,-22.484055100413200,-21.599110071782494,-20.678499539406328,-19.732580651716670,-18.769184959925067,-17.793937334878090,-16.810731632370356];
results(1).Preq15dB.optimized = [-25.499361997296640,-24.885819332865390,-24.185925781691484,-23.411186570328784,-22.571389531871922,-21.682118454936273,-20.755140676647066,-19.804694567136643];
% results(1).Preq10dB.eq_spaced = 
% results(1).Preq10dB.optimized = 

load results_8PAM
results(2).GsoadB_opt = GsoadB_opt;
results(2).Preq = Preq;
results(2).Preq20dB.eq_spaced = [-18.187839963444976,-17.199377172928640,-16.205759586855994,-15.210573129430298,-14.214582311656628,-13.218019358272807,-12.219867418954713,-11.221067726945620];
results(2).Preq20dB.optimized = [-21.462714791139070,-20.489965841775820,-19.500187454529843,-18.511590573622726,-17.526135425744190,-16.530094135385305,-15.532655576310912,-14.534177321717749];
results(2).Preq15dB.eq_spaced = [-18.113594717291104,-17.176345171343030,-16.218267971686814,-15.247568276581525,-14.268283682410010,-13.283475613259258,-12.293492197143113,-11.300610603468812];
results(2).Preq15dB.optimized = [-21.040270889040293,-20.227484272806887,-19.352627682927070,-18.445543068690355,-17.513793234886430,-16.553326162975460,-15.579540867114153,-14.596404843726539];
% results(2).Preq10dB.eq_spaced = 
% results(2).Preq10dB.optimized = 

%% Plot
legends = {};
for k = 1:length(Fn)
    figure(1), hold on, grid on, box on
    hplot(k) = plot(tx.PtxdBm, log10(ber(k).eq_spaced.count), '-o');
    legends = [legends, sprintf('Fn = %.1f dB', Fn(k))];
    
    figure(2), hold on, grid on, box on
    plot(tx.PtxdBm, log10(ber(k).optimized.count), '-o', 'Color', get(hplot(k), 'Color'));
end

for k = 1:length(Fn)
    figure(1)
    plot(tx.PtxdBm, log10(ber(k).eq_spaced.est), '-', 'Color', get(hplot(k), 'Color'));
    plot(tx.PtxdBm, log10(ber(k).eq_spaced.awgn), '--', 'Color', get(hplot(k), 'Color'));
      
    figure(2)
    plot(tx.PtxdBm, log10(ber(k).optimized.est), '-', 'Color', get(hplot(k), 'Color'));
    plot(tx.PtxdBm, log10(ber(k).optimized.awgn), '--', 'Color', get(hplot(k), 'Color'));
end

figure(1)
xlabel('Received Power (dBm)')
ylabel('log(BER)')
legend(legends{:})
axis([tx.PtxdBm([1 end]) -8 0])

figure(2)
xlabel('Received Power (dBm)')
ylabel('log(BER)')
legend(legends{:})
axis([tx.PtxdBm([1 end]) -8 0])

for k = 1:2
    figure(3), box on, grid on, hold on
    hplot = plot(Fn, results(k).GsoadB_opt.eq_spaced, '-o');
    plot(Fn, results(k).GsoadB_opt.optimized, '-s', 'Color', get(hplot, 'Color'));
    xlabel('Noise Figure (dB)')
    ylabel('Minimum SOA Gain (dB)')

    figure(4), box on, grid on, hold on
    plot(Fn, results(k).Preq.eq_spaced + results(k).GsoadB_opt.eq_spaced, '-o', 'Color', get(hplot, 'Color'));
    plot(Fn, results(k).Preq.optimized + results(k).GsoadB_opt.optimized, '-s', 'Color', get(hplot, 'Color'));
    xlabel('Noise Figure (dB)')
    ylabel('Amplifier Output Power (dBm)')
end

for k = 1:2
    figure, box on, grid on, hold on
    plot(Fn, results(k).Preq.eq_spaced, '-o');
    plot(Fn, results(k).Preq20dB.eq_spaced, '-o');
    plot(Fn, results(k).Preq15dB.eq_spaced, '-o');

    plot(Fn, results(k).Preq.optimized, '-s');
    plot(Fn, results(k).Preq20dB.optimized, '-s');
    plot(Fn, results(k).Preq15dB.optimized, '-s');
    xlabel('Noise Figure (dB)')
    ylabel('Required Transmitted Power (dBm)')
    legend('Optimal Gain', '20 dB', '15 dB')
end

