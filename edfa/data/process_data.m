%% Process data
clear, clc, close all

%% Ge:silicate
load('abs_fig1.mat')
load('gain_fig1.mat')

lamb_abs = abs_fig1(:, 1);
absorption = abs_fig1(:, 2);

lamb_gain = gain_fig1(:, 1);
gain = gain_fig1(:, 2);

clear abs_fig1 gain_fig1

% Set up fittype and options.
ft_abs = fittype( 'gauss6' );
opts_abs = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts_abs.Display = 'Off';
opts_abs.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts_abs.StartPoint = [10.0132190056334 1533.07825069569 1.71719400972773 8.48276268962462 1535.88623898282 2.19858142629229 8.395409906744 1529.90665895289 2.3017266801609 6.54849325686214 1525.53533740785 3.2752080633848 5.15236874853304 1540.68386537946 4.97526788365629 4.99014191266384 1518.69086595799 5.97939883848725];

% Fit model to data.
[fit_abs, gof] = fit( lamb_abs, absorption, ft_abs, opts_abs);

ft_gain = fittype( 'gauss6' );
opts_gain = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts_gain.Display = 'Off';
opts_gain.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts_gain.StartPoint = [9.24925662572721 1533.53846153846 1.58654543141708 7.6796121322838 1530.15384615385 2.2499968962882 7.62587887611547 1536.61538461538 2.01748525860672 5.30756302521008 1552.30769230769 3.94353897106974 5.09559145014209 1540.92307692308 3.29090439076218 4.73983202375229 1526.15384615385 4.31108447273522];

[fit_gain, gof] = fit(lamb_gain, gain, ft_gain, opts_gain);

% Plot fit with data.
figure, hold on, box on
plot(lamb_abs, absorption, 'ob')
plot(lamb_abs, fit_abs(lamb_abs), '-b')
plot(lamb_gain, gain, 'or')
plot(lamb_gain, fit_gain(lamb_gain), '-r')

save('giles_ge_silicate')

%% Al:Ge:silicate
clear
load('abs_fig2.mat')
load('gain_fig2.mat')

lamb_abs = abs_fig2(:, 1);
absorption = abs_fig2(:, 2);

lamb_gain = gain_fig2(:, 1);
gain = gain_fig2(:, 2);

clear abs_fig2 gain_fig2

% Set up fittype and options.
ft_abs = fittype( 'gauss8' );
opts_abs = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts_abs.Display = 'Off';
opts_abs.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts_abs.StartPoint = [3.77416709608372 1533.07277628032 1.42138875422518 3.44979928731032 1530.37735849057 1.65381493734226 3.14671458337035 1536.03773584906 1.95310933924247 3.11214651935932 1526.73854447439 2.00195685974451 2.35343547441453 1539.67654986523 3.13442131670349 2.33087738719468 1523.3692722372 3.2654764837632 2.00068858496203 1546.54986522911 4.88259717961085 1.93104299475082 1517.43935309973 4.66396316407884];

% Fit model to data.
[fit_abs, gof] = fit( lamb_abs, absorption, ft_abs, opts_abs);

ft_gain = fittype( 'gauss7' );
opts_gain = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts_gain.Display = 'Off';
opts_gain.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts_gain.StartPoint = [3.39559834665252 1532.92911285064 1.47054355990287 3.08307047585023 1529.81358928647 1.74080203598964 2.89057559085568 1535.62138144274 1.83761204054464 2.33824772456837 1526.41589574092 2.85383794907615 2.15742259298461 1539.38331725774 3.04479477867021 2.0791463138006 1548.70434883946 3.50901010042852 1.9275570765485 1556.67354321689 4.34904197676437];

[fit_gain, gof] = fit(lamb_gain, gain, ft_gain, opts_gain);

% Plot fit with data.
figure, hold on, box on
plot(lamb_abs, absorption, 'ob')
plot(lamb_abs, fit_abs(lamb_abs), '-b')
plot(lamb_gain, gain, 'or')
plot(lamb_gain, fit_gain(lamb_gain), '-r')

save('giles_al_ge_silicate')