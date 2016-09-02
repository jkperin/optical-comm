function y = resample_interp(x, fs1, fs2, method)
%% Resample signal x from sampling frequency fs1 to sampling frequency fs2 using linear interpolation
[~, t1] = freq_time(length(x), fs1);
dt2 = 1/fs2;
t2 = 0:dt2:t1(end);
y = interp1(t1, x, t2, method);


