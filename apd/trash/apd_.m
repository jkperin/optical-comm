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
doubly stochastic
function [Irx, v, g] = apd_doubly_stochastic(Pt, dt, tx, rx)

keff = rx.k;
Gapd = rx.Gapd;

% Constants
h = 6.62606957e-34;
q = 1.60217657e-19;
c = 299792458;

%% 1. Poisson counting process
l = Pt*dt/(h*(c/tx.lamb));  % number of photons in a interval dt with average power Pt

% quantum efficiency
eta = rx.R*(h*(c/tx.lamb))/q;

v = eta*poissrnd(l);

%% 2. Avalance gain
% Calculate random gain for every input photon
g = zeros(size(v));
for k = 1:length(v)
    vk = v(k);
    
    counter = 0;
    rmin = round(0.8*Gapd*vk);
    rmax = round(1.2*Gapd*vk);
    A = 0;
    while A < 0.9999 && counter < 50
          
        r = rmin:rmax;

        pmf = apd_gain_distribution(vk, r, keff, Gapd);
        
        A = sum(pmf);
        
        counter = counter + 1;
        
        rmin = round(0.9*rmin);
        rmax = round(1.1*rmax);
    end
    
    cdf = cumsum(pmf);
    
    u = rand(1); % uniformally-distributed r.v.
    
    [~, ix] = min(abs(cdf-u));
    
    % sample of r i.e., number of secondary electrons generated
    g(k) = r(ix);
  
end

% Detected current
Irx = (g + v)*q/dt;

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



