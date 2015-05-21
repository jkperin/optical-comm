%% Change simulation parameters and call scripts to evaluate power penalty 
%% (pp) DC-OFDM (with full and reduced dc-bias) and ACO-OFDM.

clear, clc, close all

% %% Run different cases with quantization
% for cases = 1:6
%     simulation_case = cases;
%       
%     dc_ofdm_pp_vs_f3dB_quant
%       
%     results(cases).id = sprintf('%s, CS = %d', sim.type, ofdm.CS);
%     results(cases).Fnl = Fnl;
%     results(cases).pp_ook = power_pen_ook_m;
%     results(cases).pp_awgn_ofdm = PawgndBm - PookdBm;
%     results(cases).rcliptx = RCLIPTX;
%     results(cases).rcliprx = RCLIPRX;
%     results(cases).bercount = bercount;
%     results(cases).berest = berest;  
%         
%     clear ofdm tx rx sim
%     
% end
%     
% save dc_ofdm_pp_quant results
% 
% clear results

%% Run different cases with quantization
for cases = 1:3
    simulation_case = cases;
    
    aco_ofdm_pp_vs_f3dB_quant
    
    results(cases).id = sprintf('%s, CS = %d', sim.type, ofdm.CS);
    results(cases).Fnl = Fnl;
    results(cases).pp_ook = power_pen_ook_m;
    results(cases).pp_awgn_ofdm = PawgndBm - PookdBm;
    results(cases).rcliptx = RCLIPTX;
    results(cases).rcliprx = RCLIPRX;
    results(cases).bercount = bercount;
    results(cases).berest = berest;  
    
    clear ofdm tx rx sim
    
end
   
save aco_ofdm_pp_quant results

