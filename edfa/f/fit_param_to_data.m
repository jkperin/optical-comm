function fit_param_to_data()
% Fit simulation parameter to experimental data

addpath ../
addpath ../data/
addpath ../f/
addpath ../../f/

%% Import experimental data from file
data_set_name = '-10.5dBm input per ch';

Experiment = load_experimental_data(data_set_name);

data = Experiment(data_set_name); % select desired experiment

E = EDF(6, 'corning high na');

%% 
options = optimset('Display', 'iter', 'PlotFcns', @optimplotfval);

[opt_excess_loss, ~, exitflag1] = fminbnd(@(x) iterate_excess_loss(x, E, data), 0, 1, options)

E.excess_loss = opt_excess_loss;

[opt_nsp_correction, ~, exitflag2] = fminbnd(@(x) iterate_nsp_correction(x, E, data), 1, 3, options)

    function metric = iterate_excess_loss(excess_loss, E, data)
        E.excess_loss = excess_loss;

        ReferenceBandwidth = 12.5e9; % Bandwidth for noise power measurement

        wavelength = data(1).wavelengthnm.';
        input_power_mW = dBm2Watt(data(1).InputPowerdBm.');

        Signal = Channels(wavelength * 1e-9, input_power_mW, 'forward');

        metric = 0;
        for k = 1:length(data)
            pump_power_mW = data(k).PumpPowermW;

            Pump = Channels(980e-9, pump_power_mW * 1e-3, 'forward');
            ASEf = Channels(Signal.wavelength, 0, 'forward');
            ASEb = Channels(Signal.wavelength, 0, 'backward');  

            [GaindB, Pump_out, Psignal_out, Pase, sol] = E.propagate(Pump, Signal,...
                ASEf, ASEb, ReferenceBandwidth, 'two-level', 50, false);

            PasedBm = Watt2dBm(Pase);
            
            experimental_Gain = 10.^(data(k).GaindB.'/10);
            unbiased_experimental_GaindB = 10*log10(experimental_Gain - dBm2Watt(data(k).ASEdBm.')./Signal.P);

            metric = metric + mean(abs(unbiased_experimental_GaindB - GaindB).^2)...
                            + mean(abs(data(k).ASEdBm.' - PasedBm).^2);    
        end
    end

    function metric = iterate_nsp_correction(nsp_correction, E, data)
        ReferenceBandwidth = 12.5e9; % Bandwidth for noise power measurement

        wavelength = data(1).wavelengthnm.';
        input_power_mW = dBm2Watt(data(1).InputPowerdBm.');

        Signal = Channels(wavelength * 1e-9, input_power_mW, 'forward');

        metric = 0;
        for k = 1:length(data)
            pump_power_mW = data(k).PumpPowermW;

            Pump = Channels(980e-9, pump_power_mW * 1e-3, 'forward');
            
            Pase_sa = E.analytical_ASE_PSD(Pump, Signal, nsp_correction);
            Pase_sa = Pase_sa*ReferenceBandwidth;            
            
            PasedBm = Watt2dBm(Pase_sa);

            metric = metric + mean(abs(data(k).ASEdBm.' - PasedBm).^2);    
        end        
    end
end