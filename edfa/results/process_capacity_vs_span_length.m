%% 
clear, clc, close all

addpath ../
addpath ../f/
addpath ../../f/

folder = 'capacity_vs_span_length_NL_correct';
edf_type = 'corning_type1';
ChDf = 50;
pumpWavelengthnm = 980;
Ptotal = 287*50;
Lkm = 14350;
Nspans = [478   410   359   319   287   261   239   221   205   191   179   169   159   151   144]; 
spanLengthKm = Lkm./Nspans;

for p = 1:length(Nspans)
    try
        filename = sprintf('%s/capacity_vs_span_length_EDF=%s_Ptotal=%dmW_ChDf=%dGHz_L=%dkm_Nspans=%d.mat',...
            folder, edf_type, Ptotal, ChDf, Lkm, Nspans(p)  );
        disp(filename)
        S = load(filename);

        lnm = S.lamb*1e9;

        %% Linear regime
        % Optimized power loading
        figure(201), hold on, box on
        plot(lnm, S.lin.S.PdBm)
        xlabel('Wavelength (nm)')
        ylabel('Optmized channel power (dBm)')
        if S.lin.exitflag ~= 1
            Nspans(p)
            warning('Linear regime PSO did not converge')
        end

        figure(203), hold on, box on
        hplot = plot(lnm, S.lin.num.SE);
        plot(lnm, S.lin.approx.SE, '--', 'Color', get(hplot, 'Color'))
        xlabel('Wavelength (nm)')
        ylabel('Spectral efficiency (bit/s/Hz)')
        legend('Numerical', 'Approximated')
        title('Linear regime')

        lin.SEnum(p) = sum(S.lin.num.SE);
        lin.SEapprox(p) = sum(S.lin.approx.SE);
        lin.Lopt(p) = S.lin.E.L;          

        %% Nonlinear regime
        % Optimized power loading
        figure(202), hold on, box on
        hplot = plot(lnm, S.nlin.S.PdBm, ':');
        plot(lnm, S.nlin_unc.S.PdBm, '--', 'Color', get(hplot, 'Color'))
        plot(lnm, S.nlin_sfn.S.PdBm, '-', 'Color', get(hplot, 'Color'))
        legend('PSO', 'Quasi-Newton', 'Saddle-free Newton')
        xlabel('Wavelength (nm)')
        ylabel('Optmized channel power (dBm)')

        if S.nlin.exitflag ~= 1
            Nspans(p)
            warning('Nonlinear regime PSO did not converge')
        end

        if S.nlin_unc.exitflag ~= 1
            Nspans(p)
            fprintf('Quasi-Newton exitflag = %d\n', S.nlin_unc.exitflag)
            warning('Nonlinear regime hybrid optimization did not converge')
        end
        
        if S.nlin_sfn.exitflag ~= 1
            Nspans(p)
            fprintf('Saddle-free Newton exitflag = %d\n', S.nlin_sfn.exitflag)
            warning('Nonlinear regime hybrid optimization did not converge')
        end
        
        figure(204), hold on, box on
%         hplot = plot(lnm, S.nlin.num.SE);
%         plot(lnm, S.nlin.approx.SE, '--', 'Color', get(hplot, 'Color'))
        hplot = plot(lnm, S.nlin_sfn.num.SE);
        plot(lnm, S.nlin_sfn.approx.SE, '--', 'Color', get(hplot, 'Color'))
        xlabel('Wavelength (nm)')
        ylabel('Spectral efficiency (bit/s/Hz)')
        legend('Numerical', 'Approximated')
        title('Nonlinear regime')

        nlin.SEnum(p) = sum(S.nlin.num.SE);
        nlin.SEapprox(p) = sum(S.nlin.approx.SE);
        nlin.Lopt(p) = S.nlin.E.L;

        nlin_unc.SEnum(p) = sum(S.nlin_unc.num.SE);
        nlin_unc.SEapprox(p) = sum(S.nlin_unc.approx.SE);
        nlin_unc.Lopt(p) = S.nlin_unc.E.L;
        
        nlin_sfn.SEnum(p) = sum(S.nlin_sfn.num.SE);
        nlin_sfn.SEapprox(p) = sum(S.nlin_sfn.approx.SE);
        nlin_sfn.Lopt(p) = S.nlin_sfn.E.L;
            
        drawnow

    catch e
        warning(e.message)
    end
end

figure(1), hold on, box on
plot(spanLengthKm, lin.Lopt, 'LineWidth', 2, 'DisplayName', 'Linear regime')
plot(spanLengthKm, nlin.Lopt, 'LineWidth', 2, 'DisplayName', 'Nonlinear regime')
plot(spanLengthKm, nlin_unc.Lopt, 'LineWidth', 2, 'DisplayName', 'Nonlinear regime: PSO + Quasi-Newton')
plot(spanLengthKm, nlin_sfn.Lopt, 'LineWidth', 2, 'DisplayName', 'Nonlinear regime: PSO + SFN')
xlabel('Span length (km)', 'FontSize', 12)
ylabel('Optimal EDF length (m)', 'FontSize', 12)
legend('-dynamiclegend')
set(gca, 'FontSize', 12)

figure(2), hold on, box on
hplot(1) = plot(spanLengthKm, S.problem.df*lin.SEnum/1e12, 'LineWidth', 2, 'DisplayName', 'Linear regime');
hplot(2) = plot(spanLengthKm, S.problem.df*nlin.SEnum/1e12, 'LineWidth', 2, 'DisplayName', 'Noninear regime');
hplot(3) = plot(spanLengthKm, S.problem.df*nlin_unc.SEnum/1e12, 'LineWidth', 2, 'DisplayName', 'Nonlinear regime: PSO + Quasi-Newton');
hplot(4) = plot(spanLengthKm, S.problem.df*nlin_sfn.SEnum/1e12, 'LineWidth', 2, 'DisplayName', 'Nonlinear regime: PSO + SFN');
plot(spanLengthKm, S.problem.df*lin.SEapprox/1e12, '--', 'Color', get(hplot(1), 'Color'), 'LineWidth', 2);
plot(spanLengthKm, S.problem.df*nlin.SEapprox/1e12, '--', 'Color', get(hplot(2), 'Color'), 'LineWidth', 2);
plot(spanLengthKm, S.problem.df*nlin_unc.SEapprox/1e12, '--', 'Color', get(hplot(3), 'Color'), 'LineWidth', 2);
plot(spanLengthKm, S.problem.df*nlin_sfn.SEapprox/1e12, '--', 'Color', get(hplot(4), 'Color'), 'LineWidth', 2);
legend('Linear regime', 'Nonlinear regime', 'Nonlinear regime: PSO + Quasi-Newton', 'Nonlinear regime: PSO + SFN')
xlabel('Span length (km)', 'FontSize', 12)
ylabel('Capacity per fiber (Tb/s)', 'FontSize', 12)
set(gca, 'FontSize', 12)

m = matlab2tikz(gca);
m.write_tables('capacity_vs_span_length', 'same x');