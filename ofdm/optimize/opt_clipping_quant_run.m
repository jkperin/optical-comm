clear, clc, close all

%% Main paramters

FNL = (10:5:50)*1e9;

type = {'preemphasis', 'preemphasis', 'palloc', 'palloc'};                % Type of compensation at the transmitter {'preemphasis', 'palloc'}
CS = {16 64 16 64};                       % Constellation size (effective CS in case of variable bit loading)  

%% ACO-OFDM
% Run different cases with quantization 

% for cases = [1 3]
%     sim.ENOB = 5; 
%     
%     aco_ofdm_ber_vs_clipping_ratio_rx
%         
%     clear ofdm tx rx sim ber* r*
% end
% % 
% for cases = 1:4
%     sim.ENOB = 6; 
%     
%     aco_ofdm_ber_vs_clipping_ratio_rx
%         
%     clear ofdm tx rx sim ber* r*
% end

%% DC-OFDM
% Run different cases without quantization 
for cases = [1 3]
    sim.ENOB = 5; 
    
    dc_ofdm_ber_vs_clipping_ratio_rx
        
    clear ofdm tx rx sim ber* r*
end

for cases = 1:4
    sim.ENOB = 6; 
    
    dc_ofdm_ber_vs_clipping_ratio_rx
        
    clear ofdm tx rx sim ber* r*
end