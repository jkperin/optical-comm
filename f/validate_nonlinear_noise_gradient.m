function validate_nonlinear_noise_gradient()
%% Validate gradient calculations
clear, close all

S = load('../edfa/results/capacity_vs_pump_power_PdBm/capacity_vs_pump_power_EDF=principles_type3_pump=60mW_980nm_L=286_x_50km.mat');

E = S.Eopt{S.kopt};
Signal = S.Sopt{S.kopt};
Pump = S.Pump;
problem = S.problem;
problem.step_approx = @(x) 0.5*(tanh(2*x) + 1); % Smoothing factor = 2
problem.excess_noise = 1.5; % 1.2 for 980nm, 1.6 for 1480nm
problem.epsilon = 0.05;

S = load('GN_model_coeff_spanLengthkm=50.mat');
problem.nonlinear_coeff = S.nonlinear_coeff;

P = Signal.P(60:80);
D{1} = problem.nonlinear_coeff{1}(130:170, 130:170);
D{2} = problem.nonlinear_coeff{2}(130:170, 130:170);
D{3} = problem.nonlinear_coeff{3}(130:170, 130:170);

options = optimoptions('fmincon', 'Display', 'iter', 'UseParallel', false,...
    'CheckGradients', true, 'SpecifyObjectiveGradient', true);

[PdBm, val, exitflag] = fmincon(@(PdBm) objective(PdBm, D), ...
    10*log10(P/1e-3), [], [], [], [], -20 + zeros(size(P)), zeros(size(P)), [], options);


end

function [y, dy] = objective(PdBm, D)
    % y = sum(sG.*log2(1 + SNR))
    P = 1e-3*10.^(PdBm/10);
    A = 3.801009038124394e-06;
    sG = 1:length(P); % attributes different weights to different channels
    
    [dNL, NL] = GN_model_noise_gradient(P, D);
    
    NL = NL*5^(1 + 0.05); % scaling nonlinear noise
    dNL = dNL*5^(1 + 0.05); % gradient must be scaled too
    
    SNR = P./(A + NL);
    
    y = sum(sG.*log2(1 + SNR));
    
    dSNR = diag(SNR./P) - dNL.*(SNR.^2./P);
    
    dy = 1/log(2)*sum(dSNR.*(sG./(1 + SNR)), 2);
    
    dy = log(10)/10*P.'.*dy;
end

% function [y, dy] = objective(PdBm, D)
%     % y = sum(log2(1 + SNR))
%     P = 1e-3*10.^(PdBm/10);
%     A = [3.801009038124394e-06,3.782925602748928e-06,3.740522043944979e-06,3.685112429617619e-06,3.593686001501283e-06,3.468144282694267e-06];
%     
%     [dNL, NL] = GN_model_noise_gradient(P, D);
%     
%     SNR = P./(A + NL);
%     
%     y = sum(log2(1 + SNR));
%     
%     dSNR = (SNR./P).' - dNL*(SNR.^2./P).'; % derivative of the sum of SNRs
%     dy = 1/log(2)*(dSNR./(1 + SNR.'));
%     
%     dy = log(10)/10*P.'.*dy;
% end

% function [y, dy] = objective(PdBm, D)
%     % y = sum(SNR)
%     P = 1e-3*10.^(PdBm/10);
%     A = [3.801009038124394e-06,3.782925602748928e-06,3.740522043944979e-06,3.685112429617619e-06,3.593686001501283e-06,3.468144282694267e-06];
%     
%     [dNL, NL] = GN_model_noise_gradient(P, D);
%     
%     SNR = P./(A + NL);
%     
%     y = sum(SNR);
%     
%     dy = (SNR./P).' - dNL*(SNR.^2./P).'; % derivative of the sum of SNRs
%     
%     dy = log(10)/10*P.'.*dy;
% end

% function [y, dy] = objective(PdBm, D)
%     % y = sum(NL)
%     P = 1e-3*10.^(PdBm/10);
%         
%     [dNL, NL] = GN_model_noise_gradient(P, D);
%     
%     y = sum(NL);
%     
%     dy = sum(dNL, 2); % when objective = sum(NL);   
%     dy = log(10)/10*P.'.*dy;
% end

