%% Process data from clipping ratio optimization without quantization
clear, clc, close all

ofdm_type = 'dc_ofdm';

filepath = ['results/' ofdm_type '/'];

type = {'preemphasis', 'preemphasis', 'palloc', 'palloc'};                % Type of compensation at the transmitter {'preemphasis', 'palloc'}

CS = {16 64 16 64};                       % Constellation size (effective CS in case of variable bit loading)  

FF = (10:5:50)*1e9;

rcliptx_opt = zeros(4, length(FF));

for cases = 1:4
    fprintf('----- case %d: %s, CS = %d -----\n', cases, type{cases}, CS{cases});
       
    for mm = 1:length(FF)
        sim.type = type{cases};
    
        ofdm.CS = CS{cases};
    
        filename = ['clip_opt_' ofdm_type '_' sim.type '_CS=' num2str(ofdm.CS) '_Fnl=' num2str(FF(mm)/1e9) 'GHz.mat'];
        
        load([filepath filename])
        
        rcliptx_opt(cases, mm) = ropt;
        
        clear ofdm tx rx sim ber* ropt
    end
    
    [FF.'/1e9 rcliptx_opt(cases, :).']
    
end
        
    
    
    