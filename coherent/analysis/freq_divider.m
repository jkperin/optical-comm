%% Frequency divider analysis
clear, clc, close all

addpath ../../f/
addpath ../analog/

Nsymb = 2^11;
Mct = 15;
N = Nsymb*Mct;
Rs = 56e9;
fs = Rs*Mct;
ts = 0:1/fs:(N-1)*1/fs;

phase_noise = true;
awgn_noise = true;
frequency_offset = true;

Filt = ClassFilter('bessel', 5, 0.5e9/fs, fs);
Filt.memForward = ones(size(Filt.memForward));
%
foff = 1e9;
x = cos(2*pi*foff*ts);
y = ones(size(x));

for t = Mct+1:length(x)
    [x(t), y(t-1)];
    if false % t < 100
        yt(t) = x(t);
        Filt.filter(x(t)*y(t-1));
    else
        [x(t), y(t-1)];
        y(t) = 10*Filt.filter(x(t)*y(t-1));
        1;
    end

end

figure, hold on, box on
plot(ts, x, ts, y)

% H = tf(numz, denz, 1/fs);
% woff = 2*pi*foff;
% t = 0:1/fs:2^10/fs;
% phin = pi/5;
% 
% u = sin(woff*t + phin);
% yy = lsim(H, u, t);
% 
% ytheory = (wn/woff)^2*(-sin(woff*t+phin) + sin(phin)+woff*t*cos(phin))...
%     -2*csi*wn/woff*(cos(woff*t+phin)-cos(phin));
% 
% figure
% plot(t, u, t, yy, t, ytheory, '--')

% wvec = woff;
% for k = 2:20;
%     wo = wvec(k-1) - 
%     wvec(k) = wvec(k-1) + wn^2/wvec(k-1);
% end
% 
% figure, plot(wvec/(2*pi))
    
    
    


