%% Setup edfa

addpath data/
addpath ../f

try 
    disp('1. Compile GN_model_noise.cpp and GN_model_noise_gradient.cpp')
    mex ../f/GN_model_noise.cpp
    mex ../f/GN_model_noise_gradient.cpp
    
    movefile('GN_model_noise*', '../f/')
    
catch e
    fprintf('Matlab returned error: %s.\n\n Possibly compiler is not setup\n', e.message)
    return
end

disp('2. Save EDF data as mat files')
data_corning_high_na
data_corning_new
data_corning_type1

% move save files to data/
movefile('data_corning_high_na.mat', 'data/')
movefile('data_corning_new.mat', 'data/')
movefile('data_corning_type1.mat', 'data/')

