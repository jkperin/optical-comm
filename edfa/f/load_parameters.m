function E = load_parameters(E)
%% Updates instance of EDF class with parameters corresponding to E.type
% Corning data is more complete and consistent. Other data sets were
% obtained from textbooks or papers.
switch(lower(E.type))
    %% Data provided by corning
    case 'corning_type1' 
        try
            param = load('data_corning_type1.mat');
        catch
            error('File data_corning_type1.mat not found. Add folder data/ to Matlab path.')
        end  
            
        E.excess_loss = 0; % (dB/m) excess loss due to splices for instance. 
        E.gp_980nm = 0; % gain coefficient is assumed 0, so that two-level system can be used
        E.alphap_980nm = param.pump_absorption_coeff_fun(980);% absorption cross section near 980 nm (dB/m). 
        E.core_radius = param.core_radius; % Fiber core radius
        E.doping_radius =  param.Er_radius; % Er3+ core radius
        E.rho0 = param.Er_ion_ensity; % Er3+ concentraction (cm^-3)
        E.NA = param.NA; % numerical aperture, e.g., 0.28 in [4, pg. 156]
        E.tau = param.metastable_lifetime; % metastable lifetime in s        
        E.param = param;

    case 'corning_exp'
        try
            param = load('data_corning_new.mat');
        catch
            error('File data_corning_new.mat not found. Add folder data/ to Matlab path.')
        end  
               
        E.excess_loss = 0; % (dB/m) excess loss due to splices for instance. 
        E.gp_980nm = 0; % gain coefficient is assumed 0, so that two-level system can be used
        E.alphap_980nm = param.pump_absorption_coeff_fun(980);% absorption cross section near 980 nm (dB/m). 
        E.core_radius = param.core_radius; % Fiber core radius
        E.doping_radius =  param.Er_radius; % Er3+ core radius
        E.rho0 = param.Er_ion_ensity; % Er3+ concentraction (cm^-3)
        E.NA = param.NA; % numerical aperture, e.g., 0.28 in [4, pg. 156]
        E.tau = param.metastable_lifetime; % metastable lifetime in s        
        E.param = param;
        
    case 'corning high na'
        try
            param = load('data_corning_high_na.mat');
        catch
            error('File corning_high_na.mat not found. Add folder data/ to Matlab path.')
        end  
               
        E.excess_loss = 0; % (dB/m) excess loss due to splices for instance. 
        E.gp_980nm = 0; % gain coefficient is assumed 0, so that two-level system can be used
        E.alphap_980nm = param.pump_absorption_coeff_fun(980);% absorption cross section near 980 nm (dB/m). 
        E.core_radius = param.core_radius; % Fiber core radius
        E.doping_radius =  param.Er_radius; % Er3+ core radius
        E.rho0 = param.Er_ion_ensity; % Er3+ concentraction (cm^-3)
        E.NA = param.NA; % numerical aperture, e.g., 0.28 in [4, pg. 156]
        E.tau = param.metastable_lifetime; % metastable lifetime in s        
        E.param = param;


    %% Experimental absorption and emission cross section line shapes 
    % of 1.55um transition of Er3+ in alumino-germanosilicate glass.
    % These curves correspond to Figs. 4.20, 4.21, and 4.22 of
    % "Principles and applications of EDFAs" by Desuvire
    % The crosssection curves were fit with a sum of 8 Gaussians. The 
    % coefficients below are for the peak wavelength (l), peak
    % intensity (a), and FWHM (d).
    % Fiber type I
    case 'principles_type1'
        % Fiber type I        
        edf_param.abs_peak = 7.9e-25; % Table 4.1 of Principles
        edf_param.abs_cross_sec.l = 1e-6*[1.45 1.475 1.501 1.525 1.5337 1.55 1.565 1.585];
        edf_param.abs_cross_sec.a = edf_param.abs_peak*[0.03 0.12 0.2 0.38 0.63 0.48 0.11 0.06];
        edf_param.abs_cross_sec.d = 1e-9*[30 50 40 28 10 18 20 30];      

        edf_param.ems_peak = 6.7e-25; % Table 4.1 of Principles
        edf_param.ems_cross_sec.l = 1e-6*[1.501 1.517 1.5275 1.5352 1.5515 1.5515 1.569 1.582];
        edf_param.ems_cross_sec.a = edf_param.ems_peak*[0.04 0.09 0.16 0.76 0.43 0.06 0.06 0.09];
        edf_param.ems_cross_sec.d = 1e-9*[75 40 19 8.5 20 8.5 20 50]; 

        edf_param.eta_peak = 0.84; % cross-section ratio
        
        E.param = edf_param;
        
    case 'principles_type2'    

        % Fiber type II       
        edf_param.abs_peak = 5.8e-25; % Table 4.1 of Principles
        edf_param.abs_cross_sec.l = 1e-6*[1.462 1.482 1.499 1.522 1.532 1.547 1.562];
        edf_param.abs_cross_sec.a = edf_param.abs_peak*[0.05 0.2 0.24 0.43 0.42 0.4 0.15];
        edf_param.abs_cross_sec.d = 1e-9*[30 35 36 34 13 30 30]; 

        edf_param.ems_peak = 6.5e-25; % Table 4.1 of Principles
        edf_param.ems_cross_sec.l = 1e-6*[1.462 1.499 1.525 1.533 1.553 1.561 1.586];
        edf_param.ems_cross_sec.a = edf_param.ems_peak*[0.03 0.1 0.26 0.57 0.555 0.02 0.17];
        edf_param.ems_cross_sec.d = 1e-9*[50 42 32 11 31 10 60];

        edf_param.eta_peak = 0.9; % cross-section ratio
        
        E.param = edf_param;
    case 'principles_type3' 
        % Fiber type III       
        edf_param.abs_peak = 4.7e-25; % Table 4.1 of Principles
        edf_param.abs_cross_sec.l = 1e-6*[1.44 1.482 1.492 1.515 1.53 1.544 1.555 1.57];
        edf_param.abs_cross_sec.a = edf_param.abs_peak*[0.03 0.31 0.17 0.37 0.74 0.28 0.3 0.07];
        edf_param.abs_cross_sec.d = 1e-9*[40 50 29 29 16.5 17 25 35];

        edf_param.ems_peak = 4.4e-25; % Table 4.1 of Principles
        edf_param.ems_cross_sec.l = 1e-6*[1.47 1.5 1.52 1.53 1.5425 1.556 1.575 1.6];
        edf_param.ems_cross_sec.a = edf_param.ems_peak*[0.06 0.16 0.3 0.73 0.38 0.49 0.2 0.06];
        edf_param.ems_cross_sec.d = 1e-9*[50 40 25 12.5 13 22 25 60];

        edf_param.eta_peak = 0.9;
        
        E.param = edf_param;
        
    %% Gain and absorption coefficients extracted from "Modeling Erbium-doped fiber amplifiers" by Giles and Desurvire
    % Ge:silicate fiber. Fig. 2a
    case 'giles_ge:silicate'
        try
            load('giles_ge_silicate.mat')
        catch
            error('File giles_ge_silicate.mat not found. Add folder data/ to Matlab path.')
        end
        edf_param.absorption_coeff_fun = fit_abs;
        edf_param.gain_coeff_fun = fit_gain;
        
        E.param = edf_param;
    
    % Al:Ge:silicate fiber. Fig. 2b
    case 'giles_al:ge:silicate' 
        try
            load('giles_al_ge_silicate.mat')
        catch
            error('File giles_al_ge_silicate.mat not found. Add folder data/ to Matlab path.')
        end
        edf_param.absorption_coeff_fun = fit_abs;
        edf_param.gain_coeff_fun = fit_gain;
        
        E.param = edf_param;
        
    otherwise
        error('EDF: invalid fiber type. Options are giles_ge:silicate, giles_al:ge:silicate, principles_type1, principles_type2 and principles_type3')
end