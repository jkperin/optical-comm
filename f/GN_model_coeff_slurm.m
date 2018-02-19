function D = GN_model_coeff_slurm(spanLengthkm, task)
%% Compute GN model coefficients on cluster

if any(not(isnumeric([spanLengthkm task])))
    spanLength = 1e3*str2double(spanLengthkm);
    task = round(str2double(task));
end

taskList = {[33, -1], [33, 0], [33, 1], [50, -1], [55, 0], [55, 1]};

spacing = taskList{task}(1);
l_idx = taskList{task}(2);

% Large-effective area fiber
% Assuming dispersion is 20 ps/nm/km across C-band
Fiber = fiber(spanLength, @(l) 0.165*ones(size(l)), @(l) 20.4e-6*ones(size(l)));
Fiber.gamma = 0.8e-3;

if spacing == 33
    disp('Using channel spacing of 33 GHz')
    Df = 33e9;
    dlamb = df2dlamb(Df);
    lamb = 1522e-9:dlamb:1575e-9;
elseif spacing == 50
    disp('Using channel spacing of 50 GHz')
    Df = 50e9;
    dlamb = df2dlamb(Df);
    lamb = 1522e-9:dlamb:1582e-9;
end

%% Test GN model
tic
D = GN_model_coeff(lamb, Df, Fiber, l_idx);
toc

filename = sprintf('GN_model_coeff_spanLengthkm=%skm_Df=%dGHz_l=%d.mat', spanLengthkm, spacing, l_idx);
save(filename)
