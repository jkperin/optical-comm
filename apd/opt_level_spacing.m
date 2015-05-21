function [a, b, a2, b2] = opt_level_spacing(P0, mpam, apdd, sim, N0)

a = zeros(mpam.M, 1);
b = zeros(mpam.M-1, 1);
a(1) = P0;

ther_var = N0*mpam.Rs;

% Use gaussian approximation and analytical equation
rx.N0 = N0;
rx.R = apdd.R;
rx.ka = apdd.ka;
rx.Gbw = apdd.GainBW;
rx.Gapd = apdd.Gain;
rx.Fa = @(ka, G) ka*G + (1 - ka)*(2 - 1/G);
sim.rin = 'off';
sim.shot = 'on';

tx.RIN = -Inf;

addpath f
[a2, b2] = calcOptLevelSpacing(mpam, tx, rx, sim);

for k = 1:mpam.M-1
    %% Find next decision threshold
    if a(k) ~= 0

        level_pdf = apdd.levels_pdf(a(k), 1/mpam.Rs);

        I = level_pdf.I - level_pdf.mean;

        dI = abs(I(1) - I(2));

        thermal_pdf = pdf('normal', I, 0, sqrt(ther_var));

        noise_pdf = dI*conv(level_pdf.p, thermal_pdf, 'same');
        
        noise_cdf = cumtrapz(I, noise_pdf);
        
        I = I(noise_cdf <= 1 - sim.BERtarget/2);
        
        b(k) = a(k) + I(end)/apdd.Gain;      
        
    else
        b(k) = 1/apdd.Gain*sqrt(ther_var)*qfuncinv(sim.BERtarget/2);
        1;
    end
    
    %% Find next level
    [x, ~, exitflag] = fzero(@(x) tail_prob(abs(x), b(k), apdd, ther_var, 1/mpam.Rs) - sim.BERtarget/2, a2(k+1)-b2(k));

    a(k+1) = b(k) + abs(x);
    
    if exitflag ~= 1
        warning('Function did not converge to a solution')
    end
    
%     Ptail = tail_prob(a(k+1), b(k), apdd, ther_var, 1/mpam.Rs, 'left')     
end
    
figure, hold on
plot(a, zeros(size(a)), 'xr')
plot(b, zeros(size(b)), '+k')


% Normalize
b = b/a(2);
a = a/a(2);

b2 = b2/a2(2);
a2 = a2/a2(2);


    