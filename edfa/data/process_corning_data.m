%% Process data
clear, clc, close all

%% Corning EDF
corning_edf_data

%% Fit absorption coefficient in the signal band
ft = fittype( 'gauss4' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [7.084 1528.25 3.89899497186676 5.09028646517809 1521.5 5.62878906148519 4.97020819482031 1534.75 7.32752994630869 3.8176124062692 1510 10.5576004412901];

% Fit model to data.
[fit_abs, gof] = fit(lambs, alphas, ft, opts );

c = coeffvalues(fit_abs);

fun_abs = @(x) c(1)*exp(-((x-c(2))/c(3)).^2)...
    +c(4)*exp(-((x-c(5))/c(6)).^2)...
    +c(7)*exp(-((x-c(8))/c(9)).^2)...
    +c(10)*exp(-((x-c(11))/c(12)).^2);

%% Fit gain coefficient in the signal band
% Set up fittype and options.
ft = fittype( 'gauss6' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [6.533 1530 2.34240485510081 5.55009603792048 1525.75 3.04585310236845 5.52151544868287 1534.5 3.13038287518711 4.74695277868695 1542.5 4.09169982477483 4.31239531361182 1551.5 5.24685164356404 3.72939408274568 1562 8.19502713831525];

% Fit model to data.
[fit_gain, gof] = fit( lambs, gs, ft, opts );

c = coeffvalues(fit_gain);

fun_gain = @(x) c(1)*exp(-((x-c(2))/c(3)).^2)...
    +c(4)*exp(-((x-c(5))/c(6)).^2)...
    +c(7)*exp(-((x-c(8))/c(9)).^2)...
    +c(10)*exp(-((x-c(11))/c(12)).^2)...
    +c(13)*exp(-((x-c(14))/c(15)).^2)...
    +c(16)*exp(-((x-c(17))/c(18)).^2);

% Plot fit with data.
figure, hold on, box on
plot(lambs, alphas, '.b')
plot(lambs, fit_abs(lambs), '-b')
plot(lambs, gs, '.r')
plot(lambs, fit_gain(lambs), '-r')

%% Fit absorption coefficient in the pump band
% Set up fittype and options.
ft = fittype( 'gauss5' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [4.539 978.2 1.41914574326887 4.01745803277798 975.4 1.64452437022713 3.72970783784448 980.8 1.94034087485134 3.31388995382141 972.2 2.27089684283305 2.81792029865803 984.3 3.03541530731344];

% Fit model to data.
[fit_absp, gof] = fit(lambp, alphap, ft, opts );

c = coeffvalues(fit_absp);

fun_absp = @(x) c(1)*exp(-((x-c(2))/c(3)).^2)...
    +c(4)*exp(-((x-c(5))/c(6)).^2)...
    +c(7)*exp(-((x-c(8))/c(9)).^2)...
    +c(10)*exp(-((x-c(11))/c(12)).^2)...
    +c(13)*exp(-((x-c(14))/c(15)).^2);

% Gain coefficient in the pump band is zero

figure, hold on, box on
plot(lambp, alphap, '.b')
plot(lambp, fit_absp(lambp), '-b')

% Save data
save('corning_edf')