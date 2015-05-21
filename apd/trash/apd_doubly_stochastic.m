% Calculate detected current of realistic apd
% Inputs:
% Pt = power at instant t
% dt = sampling interval of Pt
% tx = struct of transmitter (tx.lamb)
% rx = struct of apd (rx.R, rx.k, rx.Gapd)

% Outputs:
% Irx = detected current corresponding to Pt
% v = number of primary electrons generated
% g = number of secondary electrons generated

function [Irx, v, g] = apd_doubly_stochastic(Pt, dt, tx, rx)

keff = rx.ka;
Gapd = rx.Gapd;

% Constants
h = 6.62606957e-34;
q = 1.60217657e-19;
c = 299792458;

%% 1. Poisson counting process
% quantum efficiency
l = rx.R*Pt*dt/q;  % average number of primary electrons in a interval dt with average power Pt

v = poissrnd(l);

%% 2. Avalance gain
% counter = 0;
rmin = 1;
rmax = 128;
% A = 0;
% while A < 0.999 && counter < 20 % don't iterate to prevent error in gamma
% function

r = rmin:rmax;

if length(r) > 2^15
    error('apd:length', 'Vector got too large')
end

pmf = apd_gain_distribution(r, keff, Gapd);

assert(~any(isnan(pmf)), 'error while calculating pmf');
% assert(counter < 20, 'pmf calculation did not converge');

cdf = cumsum(pmf);

g = zeros(size(v));

u = linspace(0, 1, 2^14).'; % uniformly-distributed
dist = abs(bsxfun(@minus, u, cdf));
[~, ix] = min(dist, [], 2);
rp = r(ix); % r distributed accordingly to p
for kk = 1:max(v)
    g(v >= kk) = g(v >= kk) + rp(randi(2^14, 1, sum(v >= kk)));
end

% for kk = 1:length(v)
%     g(kk) = v(kk)*rp(randi(2^14, 1));
% end


% for kk = 1:max(v)
%     u = rand(sum(v >= kk), 1); % uniformly-distributed
%     dist = abs(bsxfun(@minus, u, cdf));
%     [~, ix] = min(dist, [], 2);
%     rp = r(ix); % r distributed accordingly to p
%     g(v >= kk) = g(v >= kk) + rp;
% end

% for kk = 1:length(v)
%     u = rand(v(kk), 1); % uniformly-distributed
%     dist = abs(bsxfun(@minus, u, cdf));
%     [~, ix] = min(dist, [], 2);
%     rp = r(ix); % r distributed accordingly to p
%     g(kk) = sum(rp);
% end  


% gain
fprintf('Average Gain: %.4f\n', mean(g)/mean(v))

% Detected current
Irx = g*q/dt;










%% Using generic expression for the probability of secondaries
% vu = unique(v);

% remove zero entry
% vu(vu == 0) = [];

% Calculate random gain for every input photon
% g = zeros(size(v));
% for k = 1:length(vu)
%     vk = vu(k);
%     
%     counter = 0;
%     rmin = round(0.8*(Gapd-1)*vk);
%     rmax = round(1.1*(Gapd-1)*vk);
%     A = 0;
%     while A < 0.999 && counter < 20
%           
%         r = rmin:rmax;
%         
%         if length(r) > 2^15
%             error('apd:length', 'Vector got too large')
%         end
% 
%         pmf = apd_gain_distribution(vk, r, keff, Gapd);
%         
%         A = sum(pmf);
%         
%         counter = counter + 1;
%         
%         rmin = round(0.9*rmin);
%         rmax = round(1.1*rmax);
%             
%     end
%     
%     assert(counter <= 20, 'pmf calculation did not converge');
%     
%     cdf = cumsum(pmf);
%     
%     u = rand(1, sum(vu(k) == v)); % uniformally-distributed r.v.
%     
%     ix = zeros(size(u));
%     for kk = 1:length(u)
%         [~, ix(kk)] = min(abs(cdf-u(kk)));
%     end
%     
%     % sample of r i.e., number of secondary electrons generated
%     g(vk == v) = r(ix);
%   
% end
% 
% % gain
% fprintf('Average Gain: %.4f\n', mean(g+v)/mean(v))
% 
% % Detected current
% Irx = (g + v)*q/dt;

%% Validate Random variables generated
% u = rand(1, 2^14);
% 
% for k = 1:length(u)
%     [~, ix] = min(abs(cdf-u(k)));
%     g(k) = r(ix);
% end
% 
% [G, x] = hist(g, 50);
% 
% G = G/trapz(x, G);
% 
% figure, hold on
% bar(x, G);
% plot(x, spline(r, pmf, x))
% 1;



