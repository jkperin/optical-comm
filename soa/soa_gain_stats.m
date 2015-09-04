%% Power distribution with Mch channels averaged over Ns symbols
clear, clc, close all

addpath ../mpam
addpath ../f % general functions
addpath f

dBm2W = @(x) 1e-3*10.^(x/10);
h = 6.62606957e-34; % Planck
q = 1.60217657e-19; % electron charge
c = 299792458;      % speed of light

%% SOA Parameters
carrier_lifetime = 2000e-12;
PsatdBm = 10;
G0 = 20; % unsaturated SOA gain

Psat = dBm2W(PsatdBm);

% Determine Gain x Input power curve
PfitdBm = linspace(-50, PsatdBm);
Pfit = dBm2W(PfitdBm);
Gfit = [G0 zeros(size(Pfit))];

for k = 2:length(Pfit)
    [Gfit(k), ~, exitflag] = fzero(@(G) Pfit(k) -  Psat/(abs(G)-1)*log(G0/abs(G)), Gfit(k-1));
       
    Gfit(k) = abs(Gfit(k));
    
    if exitflag ~= 1
        exitflag
    end
end

Gfit = Gfit(2:end);

figure, box on
plot(PfitdBm, 10*log10(Gfit))
xlabel('Input Power (dBm)')
ylabel('Amplifier Gain (dB)')

%% 
% M-PAM
% M, Rb, leve_spacing, pshape
mpam = PAM(4, 100e9, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));
mpam = mpam.adjust_levels(1e-3*10^(-20/10), -5);

Ns = floor(carrier_lifetime*mpam.Rs); % number of symbols to be averaged
Mch = 4; % Number of WDM channels

% 0 <= ki <= Mch x Ns 
% k1 + k2 + k3 + k4 = Mch x Ns
K = randi([0 Mch*Ns], [2^20 4]);

Ksum = K(:,1)+K(:,2)+K(:,3)+K(:,4);

K(Ksum ~= Mch*Ns, :) = [];

% Remove duplicates
k1 = 1;
k2 = 1;
while k1 < length(K)
    while k2 <= length(K) 
        if all(K(k1, :) == K(k2, :))
            K(k2, :) = [];
        end
        k2 = k2 + 1;
    end
    k1 = k1 + 1;
end

Ptot = (K(:, 1)*mpam.a(1) + K(:, 2)*mpam.a(2) + K(:, 3)*mpam.a(3) + K(:, 4)*mpam.a(4))/Ns;

figure, histfit(Ptot, 20, 'normal')

Gsoa = interp1(Pfit, Gfit, Ptot, 'spline');

figure, histfit(Gsoa, 20, 'normal')

% hist(Gsoa, 20)
% K = unique(K);


%% probability
N = 4;
r = 3;
Ni = [0 0 4;
      0 1 3;
      0 2 2;
      0 3 1;
      0 4 0;
      1 0 3;
      1 1 2;
      1 2 1;
      1 3 0;
      2 0 2;
      2 1 1;
      2 2 0;
      3 0 1;
      3 1 0;
      4 0 0];
  
 p = 1/r;
 k0 = 1;
 k1 = 2;
 sum(Ni(:, 1) ==  k0 & Ni(:, 2) == k1)/size(Ni, 1)
 
 double(k0 <= N-k1)*(nchoosek(N, k0)*p^k0*(1-p)^(N-k0))*(nchoosek(N, k1)*p^k1*(1-p)^(N-k1))