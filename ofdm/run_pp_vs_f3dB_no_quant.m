%% Change simulation parameters and call scripts to evaluate power penalty 
%% (pp) DC-OFDM (with full and reduced dc-bias) and ACO-OFDM.

clear, clc, close all

%% Run different cases without quantization
for cases = 1:4
    simulation_case = cases;
    
    ofdm.full_dc = true;
    
    dc_ofdm_pp_vs_f3dB_no_quant
      
    results(cases).id = sprintf('%s, CS = %d', sim.type, ofdm.CS);
    results(cases).Fnl = Fnl;
    results(cases).pp_ook = power_pen_ook_m;
    results(cases).rcliptx = RCLIP;
    results(cases).bercount = bercount;
    results(cases).berest = berest;  
    results(cases).pp_awgn_ofdm = PawgndBm - PookdBm;
    
    clear ofdm tx rx sim
    
end
    
save dc_ofdm_pp_full_dc results

clear results

%% Run different cases without quantization
for cases = 1:4
    simulation_case = cases;
    
    ofdm.full_dc = false;
    
    dc_ofdm_pp_vs_f3dB_no_quant
      
    results(cases).id = sprintf('%s, CS = %d', sim.type, ofdm.CS);
    results(cases).Fnl = Fnl;
    results(cases).pp_ook = power_pen_ook_m;
    results(cases).rcliptx = RCLIP;
    results(cases).bercount = bercount;
    results(cases).berest = berest;  
    results(cases).pp_awgn_ofdm = PawgndBm - PookdBm;
    
    clear ofdm tx rx sim
    
end
    
save dc_ofdm_pp results

clear results

%% Run different cases without quantization
for cases = 1:4
    simulation_case = cases;
    
    aco_ofdm_pp_vs_f3dB_no_quant
    
    results(cases).id = sprintf('%s, CS = %d', sim.type, ofdm.CS);
    results(cases).Fnl = Fnl;
    results(cases).pp_ook = power_pen_ook_m;
    results(cases).rcliptx = RCLIP;
    results(cases).bercount = bercount;
    results(cases).berest = berest;  
    results(cases).pp_awgn_ofdm = PawgndBm - PookdBm;
    
    clear ofdm tx rx sim
    
end
   
save aco_ofdm_pp results

