function Experiment = load_experimental_data(Sheets, enable_plot)

folder = 'data/';
filename = 'Amp2Characterization.xlsx';

% Spreadsheet details
WavelengthRange = "B3:B42";
ReferencePowerRange = "C3:C42";
PumpPowerRanges = ["F1", "J1", "N1", "R1", "V1"];
OutputPowerRanges = ["E3:E42", "I3:I42", "M3:M42", "Q3:Q42", "U3:U42"];
ASERanges = ["F3:F42", "J3:J42", "N3:N42", "R3:R42", "V3:V42"];
GainRanges = ["G3:G42", "K3:K42", "O3:O42", "S3:S42", "W3:W42"];

if not(exist('Sheets', 'var'))
    Sheets = ["-16.5dBm input per ch", "-14.5dBm input per ch", "-12.5dBm input per ch",...
        "-10.5dBm input per ch"];
else
    if class(Sheets) == 'char'
        Sheets = string(Sheets);
    end
end

% Signal total power was measured at point A before the splice.  
SignalAttdB = 0.05 + 0.15 + 0.25; % Splice loss at A + Splice loss at B + Loss of WDM coupler
PumpAttdB = 0.15; % Splice loss at B
GainOffsetdB = 0.35; % Amount by which gain was underestimated dB
MeasurementLoss = 0.15 + 0.7 + 1.8; % Splice loss at C + Loss of switch + Loss at input to OSA 

Experiment = containers.Map();
figure_count = 0;
for sheet = Sheets
    wavelengthnm = xlsread([folder filename], sheet, WavelengthRange);
    RefInputPowerdBm = xlsread([folder filename], sheet, ReferencePowerRange);
    
    for k = 1:length(GainRanges) 
        data(k).wavelengthnm = wavelengthnm;
        data(k).InputPowerdBm = RefInputPowerdBm;
        data(k).NominalPumpPowermW = xlsread([folder filename], sheet, PumpPowerRanges(k));
        data(k).OutputPowerdBm = xlsread([folder filename], sheet, OutputPowerRanges(k));
        data(k).ASEdBm = xlsread([folder filename], sheet, ASERanges(k));
        data(k).GaindB = xlsread([folder filename], sheet, GainRanges(k));
        
        % Adjust raw measurents
        data(k).PumpPowermW = dBm2Watt(Watt2dBm(data(k).NominalPumpPowermW) - PumpAttdB);
        data(k).InputPowerdBm = data(k).InputPowerdBm - SignalAttdB + MeasurementLoss;
        data(k).GaindB = data(k).GaindB + GainOffsetdB;
        data(k).ASEdBm = data(k).ASEdBm + MeasurementLoss;
    
        if exist('enable_plot', 'var') && enable_plot
            leg = sprintf('Ppump = %.2f mW', data(k).PumpPowermW);
            figure(80 + figure_count)
            subplot(311), hold on
            plot(data(k).wavelengthnm, data(k).InputPowerdBm, 'DisplayName', leg)
            xlabel('Wavelength (nm)')
            ylabel('Input power (dBm)')
            legend('-dynamiclegend')
            title(sheet)

            subplot(312), hold on
            plot(data(k).wavelengthnm, data(k).GaindB, 'DisplayName', leg)
            xlabel('Wavelength (nm)')
            ylabel('Gain (dB)')
            legend('-dynamiclegend')

            subplot(313), hold on
            plot(data(k).wavelengthnm, data(k).ASEdBm, 'DisplayName', leg)
            xlabel('Wavelength (nm)')
            ylabel('ASE (dBm)')
            legend('-dynamiclegend')

            drawnow
        end
    end
   
    Experiment(char(sheet)) = data;
    figure_count = figure_count + 1;
end