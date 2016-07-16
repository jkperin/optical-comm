%% BER vs OSNR for 56 Gbps 4-PAM system in back to back
clear, clc, close all

folder = 'data/waveforms/b2b_OSNR/';
files = ls([folder '*.h5']);
dacfile = 'data/waveforms/b2b_OSNR/pam4_rect_Rb=55Gbps_preemph';

for k = 1:size(files, 1)
    file = [folder files(k, :)];
    
    OSNR(k) = str2double(file(strfind(file, 'osnr=')+5:strfind(file, 'dB')-1));
    
    ber_count(k) = process_pam_waveforms(file, dacfile);
end

[OSNR, idx] = sort(OSNR);
ber_count = ber_count(idx);

figure, box on
plot(OSNR, log10(ber_count), '-ok', 'linewidth', 2)
xlabel('OSNR (dB)', 'FontSize', 12)
ylabel('log_{10}(BER)', 'FontSize', 12)
set(gca, 'FontSize', 12)
title('BER vs OSNR for 56 Gbps 4-PAM system in back to back')
axis([23 42 -6.5 -2.5])
    
    
    