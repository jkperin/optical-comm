% Combine GN model results
clear, clc, close all

FILE = 'GN_model_coeff_spanLengthkm=50_Df=33GHz';

load([FILE '_l=-1.mat'])
nonlinear_coeff{1} = D;

load([FILE '_l=0.mat']);
nonlinear_coeff{2} = D;

load([FILE '_l=1.mat']);
nonlinear_coeff{3} = D;

clear D l_idx task taskList spanLengthkm
save(FILE)