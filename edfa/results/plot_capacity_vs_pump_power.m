%% 
clear, clc, close all

addpath ../
addpath ../f/
addpath ../../f/

% paper results
% folder =  'capacity_vs_pump_power_highNA_final';
% pumpPowermW = [25:5:200 225:25:400];
% Nspans = 287; % 317

folder = 'edf_6m_1dB_margin';
pumpPowermW = [40 60 80];
Nspans = 220;

edf_type = 'corning high NA';
ChDf = 50;
pumpWavelengthnm = 980;
spanLengthKm = 50; % 

for p = 1:length(pumpPowermW)
    try
        filename = sprintf('%s/capacity_vs_pump_power_EDF=%s_pump=%dmW_%dnm_ChDf=%dGHz_L=%d_x_%dkm.mat',...
            folder, edf_type, pumpPowermW(p), pumpWavelengthnm, ChDf, Nspans, spanLengthKm);
        disp(filename)
        S = load(filename);     

        lin.SEnum(p) = sum(S.lin.num.SE);
        lin.SEapprox(p) = sum(S.lin.approx.SE);
        lin.Lopt(p) = S.lin.E.L;          

        if S.nlin.exitflag ~= 1
            pumpPowermW(p)
            warning('Nonlinear regime PSO did not converge')
        end

        if S.nlin_unc.exitflag ~= 1
            pumpPowermW(p)
            warning('Nonlinear regime hybrid optimization did not converge')
        end
        
        if S.nlin_sfn.exitflag ~= 1
            pumpPowermW(p)
            warning('Nonlinear regime SFN optimization did not converge')
        end

        nlin.SEnum(p) = sum(S.nlin.num.SE);
        nlin.SEapprox(p) = sum(S.nlin.approx.SE);
        nlin.Lopt(p) = S.nlin.E.L;

        nlin_unc.SEnum(p) = sum(S.nlin_unc.num.SE);
        nlin_unc.SEapprox(p) = sum(S.nlin_unc.approx.SE);
        nlin_unc.Lopt(p) = S.nlin_sfn.E.L;

        nlin_sfn.SEnum(p) = sum(S.nlin_sfn.num.SE);
        nlin_sfn.SEapprox(p) = sum(S.nlin_sfn.approx.SE);
        nlin_sfn.Lopt(p) = S.nlin_sfn.E.L;
        
        InputPower(p) = Watt2dBm(sum(S.nlin_sfn.S.P));
        LaunchedPower(p) = Watt2dBm(sum(S.nlin_sfn.S.P)) + S.problem.spanAttdB;
        NLpower(p) = sum(S.nlin_sfn.num.NL(S.nlin_sfn.S.P ~= 0));
        ASEpower(p) = sum(S.nlin_sfn.num.Pase(S.nlin_sfn.S.P ~= 0));
        
        SignalOut = S.lin.S;
        SignalOut.PdBm = SignalOut.PdBm + S.lin.num.GaindB;
        PCElin(p) = S.lin.E.power_conversion_efficiency(S.Pump, S.lin.S, SignalOut);  
        
        SignalOut = S.nlin_sfn.S;
        SignalOut.PdBm = SignalOut.PdBm + S.nlin_sfn.num.GaindB;
        PCEnlin(p) = S.nlin_sfn.E.power_conversion_efficiency(S.Pump, S.nlin_sfn.S, SignalOut);
        
        maxSignalP(p) = max(S.nlin_sfn.S.P);
        Pcap(p) = S.problem.Pon;
        
    catch e
        warning(e.message)
    end
end

figure(1), hold on, box on
plot(pumpPowermW, lin.Lopt, 'LineWidth', 2, 'DisplayName', 'Linear regime')
plot(pumpPowermW, nlin.Lopt, 'LineWidth', 2, 'DisplayName', 'Nonlinear regime')
plot(pumpPowermW, nlin_unc.Lopt, 'LineWidth', 2, 'DisplayName', 'Nonlinear regime: PSO + Quasi-Newton')
plot(pumpPowermW, nlin_sfn.Lopt, 'LineWidth', 2, 'DisplayName', 'Nonlinear regime: PSO + SFN')
xlabel('Pump power (mW)', 'FontSize', 12)
ylabel('Optimal EDF length (m)', 'FontSize', 12)
legend('-dynamiclegend')
set(gca, 'FontSize', 12)

figure(2), hold on, box on
hplot(1) = plot(pumpPowermW, S.problem.df*lin.SEnum/1e12, 'LineWidth', 2, 'DisplayName', 'Linear regime');
hplot(2) = plot(pumpPowermW, S.problem.df*nlin.SEnum/1e12, 'LineWidth', 2, 'DisplayName', 'Noninear regime');
hplot(3) = plot(pumpPowermW, S.problem.df*nlin_unc.SEnum/1e12, 'LineWidth', 2, 'DisplayName', 'Nonlinear regime: PSO + Quasi-Newton');
hplot(4) = plot(pumpPowermW, S.problem.df*nlin_sfn.SEnum/1e12, 'LineWidth', 2, 'DisplayName', 'Nonlinear regime: PSO + SFN');
plot(pumpPowermW, S.problem.df*lin.SEapprox/1e12, '--', 'Color', get(hplot(1), 'Color'), 'LineWidth', 2);
plot(pumpPowermW, S.problem.df*nlin.SEapprox/1e12, '--', 'Color', get(hplot(2), 'Color'), 'LineWidth', 2);
plot(pumpPowermW, S.problem.df*nlin_unc.SEapprox/1e12, '--', 'Color', get(hplot(3), 'Color'), 'LineWidth', 2);
plot(pumpPowermW, S.problem.df*nlin_sfn.SEapprox/1e12, '--', 'Color', get(hplot(4), 'Color'), 'LineWidth', 2);
legend('Linear regime', 'Nonlinear regime', 'Nonlinear regime: PSO + Quasi-Newton', 'Nonlinear regime: PSO + SFN')
xlabel('Pump power (mW)', 'FontSize', 12)
ylabel('Capacity per fiber (Tb/s)', 'FontSize', 12)
set(gca, 'FontSize', 12)
axis([20 300 10 50])
% 
m = matlab2tikz(gca);
m.write_tables('capacity_vs_pump_power', 'same x')

%% ASE vs nonlinear power
figure(3), hold on, box on
plot(pumpPowermW, Watt2dBm(ASEpower) - Watt2dBm(NLpower), 'LineWidth', 2)
xlabel('Pump power (mW)', 'FontSize', 12)
ylabel('ASE to NL noise ratio (dB)', 'FontSize', 12)
set(gca, 'FontSize', 12)

m = matlab2tikz(gca);
m.write_tables('ASE_vs_NL', 'same x')
    
%% PCE
figure(4), hold on, box on
plot(pumpPowermW, 100*PCElin, 'LineWidth', 2)
plot(pumpPowermW, 100*PCEnlin, 'LineWidth', 2)
xlabel('Pump power (mW)', 'FontSize', 12)
ylabel('Power conversion efficiency (%)', 'FontSize', 12)
set(gca, 'FontSize', 12)
legend('Linear', 'Nonlinear')
m = matlab2tikz(gca);
m.write_tables('PCE_vs_Ppump', 'same x')

%% Launched power
figure, hold on, box on
plot(pumpPowermW, InputPower)
plot(pumpPowermW, LaunchedPower)
xlabel('Pump power (mW)', 'FontSize', 12)
ylabel('Power (dBm)', 'FontSize', 12)
legend('Input', 'Launched')
set(gca, 'FontSize', 12)

m = matlab2tikz(gca);
m.write_tables('optical_power_vs_pump', 'same x')

%%
Cnum = S.problem.df*nlin_sfn.SEnum/1e12;
fo = fitoptions('Method','NonlinearLeastSquares',...
           'Lower',[0,1, -Inf -10],...
           'Upper',[Inf,Inf, Inf Inf],...
           'StartPoint',[1 1 1 0]);
ft = fittype('a*log2(d + b*x) + c','options',fo);
Cfit = fit(pumpPowermW.', Cnum.', ft);
pp = linspace(pumpPowermW(1), pumpPowermW(end));
figure, hold on
plot(pumpPowermW, S.problem.df*nlin_sfn.SEnum/1e12)
hplt = plot(pp, Cfit(pp), ':k', 'LineWidth', 2);

V = 12e3;
Nspans = 287;
Ptot = V^2/(4*1*Nspans*spanLengthKm); % total electrical power 
Po =  0:0.1:0.3; % power spent in other operations

maxDim = 40;
Ppump = zeros(length(Po), maxDim);
Cp = zeros(length(Po), maxDim);
leg = {};
for p = 1:length(Po)
    for s = 1:maxDim
        Ppump(p, s) = max(0.4*(Ptot/(2*s*Nspans)-Po(p)), 0); % optical power per edfa assuming 1 spatial dimension
        Cp(p, s) = s*Cfit(Ppump(p, s)*1e3);
    end
    leg = [leg sprintf('P_o = %.1f W', Po(p))];
    
    Cp(p, Ppump(p, :) < 20e-3) = NaN;
end

Cp(Cp < 0) = NaN;
% Cp = Cp/Cp(1, 1);
figure, hold on, box on
plot(1:maxDim, Cp)
%plot(1:maxDim, (1:maxDim)/2, ':k')
legend(leg)
set(gca, 'FontSize', 12)
xlabel('Spatial dimensions', 'FontSize', 12)
ylabel('Cable capacity (Tb/s)', 'FontSize', 12)
legend('-dynamiclegend')


[Cmax, idx] = max(Cp, [], 2)
PpumpmW = Ppump*1e3;
PpumpmW(1, idx(1))
PpumpmW(2, idx(2))
PpumpmW(3, idx(3))
PpumpmW(4, idx(4))

m = matlab2tikz(gca);
m.write_tables('capacity_vs_spatial_dims', 'same x')