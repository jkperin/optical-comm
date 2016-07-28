%% BER vs OSNR for 56 Gbps 4-PAM system in back to back
clear, clc

folder = 'data/waveforms/BER_vs_OSNR_b2b/set1/';
Rb = 56e9;
Rs = 56e9/2;
% files = {'data/waveforms/b2b_OSNR/wave-4PAM-56Gbps-preemph-0km-osnr=23dB.h5',...
% 'data/waveforms/b2b_OSNR/wave-4PAM-56Gbps-preemph-0km-osnr=27dB.h5',...
% 'data/waveforms/b2b_OSNR/wave-4PAM-56Gbps-preemph-0km-osnr=29dB.h5',...
% 'data/waveforms/b2b_OSNR/wave-4PAM-56Gbps-preemph-0km-osnr=31dB.h5',...
% 'data/waveforms/b2b_OSNR/wave-4PAM-56Gbps-preemph-0km-osnr=33dB.h5',...
% 'data/waveforms/b2b_OSNR/wave-4PAM-56Gbps-preemph-0km-osnr=35dB.h5',...
% 'data/waveforms/b2b_OSNR/wave-4PAM-56Gbps-preemph-0km-osnr=38dB.h5',...
% 'data/waveforms/b2b_OSNR/wave-4PAM-56Gbps-preemph-0km-osnr=40dB.h5',...
% 'data/waveforms/b2b_OSNR/wave-4PAM-56Gbps-preemph-0km-osnr=42dB.h5',...
% 'data/waveforms/b2b_OSNR/wave-4PAM-56Gbps-preemph-0km-osnr=44dB.h5',...
% 'data/waveforms/b2b_OSNR/wave-4PAM-56Gbps-preemph-0km-osnr=46dB.h5',...
% 'data/waveforms/b2b_OSNR/wave-4PAM-56Gbps-preemph-0km-osnr=48dB.h5'};
% dacfile = 'data/waveforms/b2b_OSNR/pam4_rect_Rb=55Gbps_preemph';

%% New experiment (set 1)
files = {'wave-4PAM-56Gbps-preemph-0km-osnr=25dB.h5',...
'wave-4PAM-56Gbps-preemph-0km-osnr=26dB.h5',...
'wave-4PAM-56Gbps-preemph-0km-osnr=27dB.h5',...
'wave-4PAM-56Gbps-preemph-0km-osnr=28dB.h5',...
'wave-4PAM-56Gbps-preemph-0km-osnr=29dB.h5',...
'wave-4PAM-56Gbps-preemph-0km-osnr=30dB.h5',...
'wave-4PAM-56Gbps-preemph-0km-osnr=31dB.h5',...
'wave-4PAM-56Gbps-preemph-0km-osnr=32dB.h5',...
'wave-4PAM-56Gbps-preemph-0km-osnr=34dB.h5',...
'wave-4PAM-56Gbps-preemph-0km-osnr=36dB.h5',...
'wave-4PAM-56Gbps-preemph-0km-osnr=39dB.h5',...
'wave-4PAM-56Gbps-preemph-0km-osnr=43dB.h5',...
'wave-4PAM-56Gbps-preemph-0km-osnr=48dB.h5'};
dacfile = 'data/waveforms/BER_vs_OSNR_b2b/pam4_rect_Rb=55Gbps_preemph.mat';

%% New experiment (set 2)
% files = {'wave-4PAM-56Gbps-preemph-0km-osnr=26dB.h5',...
% 'wave-4PAM-56Gbps-preemph-0km-osnr=28dB.h5',...
% 'wave-4PAM-56Gbps-preemph-0km-osnr=30dB.h5',...
% 'wave-4PAM-56Gbps-preemph-0km-osnr=32dB.h5',...  
% 'wave-4PAM-56Gbps-preemph-0km-osnr=34dB.h5',...  
% 'wave-4PAM-56Gbps-preemph-0km-osnr=36dB.h5',...  
% 'wave-4PAM-56Gbps-preemph-0km-osnr=38dB.h5',...  
% 'wave-4PAM-56Gbps-preemph-0km-osnr=41dB.h5',...  
% 'wave-4PAM-56Gbps-preemph-0km-osnr=46dB.h5 '};
% dacfile = 'data/waveforms/BER_vs_OSNR_b2b/pam4_rect_Rb=55Gbps_preemph.mat';

%% Predistortion with Pswing = 0.8, and driver Vg = 2V
% files = {'wave-4PAM-56Gbps-preemph_predist-0km-osnr=26dB.h5',...
%     'wave-4PAM-56Gbps-preemph_predist-0km-osnr=28dB.h5',...
%     'wave-4PAM-56Gbps-preemph_predist-0km-osnr=30dB.h5',...
%     'wave-4PAM-56Gbps-preemph_predist-0km-osnr=32dB.h5',...
%     'wave-4PAM-56Gbps-preemph_predist-0km-osnr=34dB.h5',...
%     'wave-4PAM-56Gbps-preemph_predist-0km-osnr=36dB.h5',...
%     'wave-4PAM-56Gbps-preemph_predist-0km-osnr=38dB.h5',...
%     'wave-4PAM-56Gbps-preemph_predist-0km-osnr=40dB.h5',...
%     'wave-4PAM-56Gbps-preemph_predist-0km-osnr=41dB.h5',...
%     'wave-4PAM-56Gbps-preemph_predist-0km-osnr=45dB.h5'};
% dacfile = 'pam4_rect_Rb=55Gbps_preemph_predist.mat';

parfor k = 1:length(files)  
    file = files{k}
    
    OSNRdB(k) = str2double(file(strfind(file, 'osnr=')+5:strfind(file, 'dB')-1));
    
    ber_count(k) = process_pam_waveforms([folder file], dacfile);
end

% [OSNRdB, idx] = sort(OSNRdB);
% ber_count = ber_count(idx);
% OSNRdBtheory = linspace(OSNRdB(1), OSNRdB(end));
% 
% figure(1), box on, hold on
% plot(OSNRdB, log10(ber_count), '-ok', 'linewidth', 2)
% plot(OSNRdBtheory, log10(pam_ber_from_osnr(4, OSNRdBtheory, Rs/2)), ':k', 'linewidth', 2)
% xlabel('OSNR (dB)', 'FontSize', 12)
% ylabel('log_{10}(BER)', 'FontSize', 12)
% set(gca, 'FontSize', 12)
% title('BER vs OSNR for 56 Gbps 4-PAM system in back to back')
% axis([23 42 -8 -2])
    
    
    