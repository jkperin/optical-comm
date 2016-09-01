clear, clc, close all

addpath ../../mpam/
load('symb_centers')

% Points
centers(:, 1) = real(symbol_centers_64qam(:, 1));
centers(:, 2) = imag(symbol_centers_64qam(:, 1));
centers(:, 3) = real(symbol_centers_64qam(:, 2));
centers(:, 4) = imag(symbol_centers_64qam(:, 2));

mzm_nonlinearity = @(levels, V) V(3)*sin(pi/2*(levels*V(1) + V(2))) + V(4);

M = 64;
mpam = PAM(sqrt(M), 1);
dlevel = 2/(sqrt(M)-1);
ideal_levels = -1:dlevel:1;
mpam.a = ideal_levels.';
titles = {'Pol X, I', 'Pol X, Q', 'Pol Y, I','Pol Y, Q'};

for k = 1:4
    centerk = centers(:, k).';
    dlevel = (max(centerk)-min(centerk))/(sqrt(M)-1);
    thresholds = min(centerk)+dlevel/2:dlevel:max(centerk);
    mpam.b = thresholds.';
    levels = mpam.mod(mpam.demod(centerk));
    
    [Vset, fval, exitflag] = fminsearch(@(V) norm(centerk...
        - (mzm_nonlinearity(levels, V) - mean(mzm_nonlinearity(levels, V)))), [0.5 -0.1 1 0]);

    fval
    if exitflag ~= 1
        disp('4-PAM bias control did not converge')
    end

    fprintf('====== %s ======\n', titles{k})
    fprintf('Vswing/Vpi = %.2f\n', 2*Vset(1))
    fprintf('Vbias/Vpi = %.2f\n', Vset(2))

    Vdrive = ideal_levels*Vset(1) + Vset(2);
    Pset = sin(pi/2*Vdrive);

    figure(1);
    subplot(2, 2, k), hold on, box on
    t = linspace(-1, 1);
    plot(t, sin(pi/2*t), 'k');
    plot(([1; 1]*Vdrive), [zeros(1, mpam.M)-1; Pset], 'k');
    plot([zeros(1, mpam.M)-1; Vdrive], ([1; 1]*Pset), 'k');
    xlabel('Driving signal')
    ylabel('Resulting power levels')
    axis([-1 1 -1 1])
    title(titles{k})   
    drawnow
    
    figure(2)
    subplot(2, 2, k), hold on, box on
    stem(1:mpam.M-1, diff(Pset))
    ylabel('Optical level spacing')  
    xlabel('Level ')
    set(gca, 'xtick', 1:mpam.M-1)
    title(titles{k})  
    axis([1, mpam.M-1 0 0.12])
end