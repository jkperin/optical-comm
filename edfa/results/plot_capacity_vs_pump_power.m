%% 
clear, clc, close all

addpath ../
addpath ../f/
addpath ../../f/

folder = 'capacity_vs_pump_power';
edf_type = 'corning_type1';
ChDf = 50;
pumpWavelengthnm = 980;
pumpPowermW = [20:5:150 160:10:200 250:50:400]; %%[30:5:100 150:50:250 275:25:400 450:50:1000];
% pumpPowermW = [50 150];
Nspans = 287; % 317
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
Cnum = S.problem.df*nlin_sfn.SEnum/1e12;
fo = fitoptions('Method','NonlinearLeastSquares',...
           'Lower',[0,0,-Inf],...
           'Upper',[Inf,Inf, Inf],...
           'StartPoint',[1 1 0]);
ft = fittype('a*log2(1 + b*x) + c','options',fo);
Cfit = fit(pumpPowermW.', Cnum.', ft);
hplot(1) = plot(pumpPowermW, S.problem.df*lin.SEnum/1e12, 'LineWidth', 2, 'DisplayName', 'Linear regime');
hplot(2) = plot(pumpPowermW, S.problem.df*nlin.SEnum/1e12, 'LineWidth', 2, 'DisplayName', 'Noninear regime');
hplot(3) = plot(pumpPowermW, S.problem.df*nlin_unc.SEnum/1e12, 'LineWidth', 2, 'DisplayName', 'Nonlinear regime: PSO + Quasi-Newton');
hplot(4) = plot(pumpPowermW, S.problem.df*nlin_sfn.SEnum/1e12, 'LineWidth', 2, 'DisplayName', 'Nonlinear regime: PSO + SFN');
% plot(pumpPowermW, S.problem.df*lin.SEapprox/1e12, '--', 'Color', get(hplot(1), 'Color'), 'LineWidth', 2);
% plot(pumpPowermW, S.problem.df*nlin.SEapprox/1e12, '--', 'Color', get(hplot(2), 'Color'), 'LineWidth', 2);
% plot(pumpPowermW, S.problem.df*nlin_unc.SEapprox/1e12, '--', 'Color', get(hplot(3), 'Color'), 'LineWidth', 2);
% plot(pumpPowermW, S.problem.df*nlin_sfn.SEapprox/1e12, '--', 'Color', get(hplot(4), 'Color'), 'LineWidth', 2);
hplt = plot(pumpPowermW, Cnum, ':', 'LineWidth', 2);
legend('Linear regime', 'Nonlinear regime', 'Nonlinear regime: PSO + Quasi-Newton', 'Nonlinear regime: PSO + SFN')
xlabel('Pump power (mW)', 'FontSize', 12)
ylabel('Capacity per fiber (Tb/s)', 'FontSize', 12)
set(gca, 'FontSize', 12)
% 
m = matlab2tikz(gca);
m.write_tables('capacity_vs_pump_power', 'same x')

V = 12e3;
Nspans = 287;
I = V/(2*1*Nspans*spanLengthKm);
Ptot = V*I/2; % total electrical power 
Po = 0.2; % power spent in other operations

Po =  0:0.1:0.5; % power spent in other operations
leg = {};
for p = 1:length(Po)
    for s = 1:20
        Ppump = max(0.4*(Ptot/(2*s*Nspans)-Po(p)), 0); % optical power per edfa assuming 1 spatial dimension
        Cp(p, s) = s*Cfit(Ppump*1e3);
    end
    leg = [leg sprintf('P_o = %.1f W', Po(p))];
end

Cp(Cp < 0) = NaN;
Cp = Cp/Cp(1, 1);
figure, hold on, box on
plot(1:20, Cp)
legend(leg)
set(gca, 'FontSize', 12)
xlabel('Spatial dimensions', 'FontSize', 12)
ylabel('Capacity improvement', 'FontSize', 12)
legend('-dynamiclegend')

m = matlab2tikz(gca);
m.write_tables('capacity_vs_spatial_dims', 'same x')