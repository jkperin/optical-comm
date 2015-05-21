clear, clc, close all

format compact

addpath f           % functions path

rng('default')      % initiate default random number generator
rng('shuffle');     % Reinitialize the random number generator used by rand, randi, and randn with a seed based on the current time

sim.Navg = 5;

%% CS = 64
cases.ENOB = {6, 6,... % no shot no RIN
              6, 6,... % shot + RIN
              6, 6,... % only shot
              6, 6};   % only RIN

cases.type = {'preemphasis', 'palloc',...
              'preemphasis', 'palloc',...
              'preemphasis', 'palloc',...
              'preemphasis', 'palloc'};
          
cases.shot = {false, false,...
              true, true,...
              true, true,...
              false, false};

cases.RIN = {false, false,...
              true, true,...
              false, false,...
              true, true};

          
for ncase = 1:4
   
    ofdm.CS = 64;
          
    sim.type = cases.type{ncase};                % type of power allocation

    sim.ENOB = cases.ENOB{ncase};                       % Effective number of bits for DAC and ADC (only used if sim.quantiz is true)

    sim.shot = cases.shot{ncase};                    % Include shot noise?

    sim.RIN = cases.RIN{ncase};                     % Include intensity noise?
    
    dc_ofdm_ber_vs_L
    
    results_CS_64QAM{ncase}.berc = bercount;
    results_CS_64QAM{ncase}.bere = berest;
    results_CS_64QAM{ncase}.power_pen_ook_m = power_pen_ook_m;
    results_CS_64QAM{ncase}.power_pen_ook_e = power_pen_ook_e;
end

save('results_CS_64QAM', 'results_CS_64QAM', 'cases')