clear, clc, close all

%% Main paramters

FNL = (10:5:50)*1e9;

type = {'preemphasis', 'preemphasis', 'palloc', 'palloc'};                % Type of compensation at the transmitter {'preemphasis', 'palloc'}
CS = {16 64 16 64};                       % Constellation size (effective CS in case of variable bit loading)  

%% ACO-OFDM
% Run different cases without quantization 
for cases = 1:4
    
    fprintf('----- case %d: %s, CS = %d -----\n', cases, type{cases}, CS{cases});
    
    sim.type = type{cases};
    ofdm.CS = CS{cases};
    
    for ll = 1:length(FNL)
        Fnl = FNL(ll);
        
        aco_ofdm_ber_vs_clipping_ratio
    end
    
    clear ofdm tx rx sim ber* r*
end

%% DC-OFDM
% Run different cases without quantization 
for cases = 1:4
    
    fprintf('----- case %d: %s, CS = %d -----\n', cases, type{cases}, CS{cases});
    
    sim.type = type{cases};
    ofdm.CS = CS{cases};
    
    for ll = 1:length(FNL)
        Fnl = FNL(ll);
        
        dc_ofdm_ber_vs_clipping_ratio
    end
    
    clear ofdm tx rx sim ber* r*
end