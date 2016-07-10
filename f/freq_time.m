function [f, t] = freq_time(N, fs)
%% Generate frequency and time measures
dt = 1/fs;
t = 0:dt:(N-1)*dt;
df = 1/(dt*N);
f = -fs/2:df:fs/2-df;