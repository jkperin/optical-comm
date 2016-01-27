% scamp2.m   Joseph M. Kahn  12-2-11
% Models amplification of optical field by a scalar system having frequency-dependent gain.
function Eamp = scamp(Ein,deltat,G,S)

% x and y components of input signal
Einx = Ein(1,:);
Einy = Ein(2,:);

% generate independent, identically distributed noises for x and y, I and Q, each with variance  2/deltat
Nx = sqrt(1/2/deltat)*(randn(size(Einx)) + 1i*randn(size(Einx)));
Ny = sqrt(1/2/deltat)*(randn(size(Einy)) + 1i*randn(size(Einy)));

% scale signal components to amplify by frequency-dependent gain G
Eampsx = ifft(sqrt(G).*fft(Einx));
Eampsy = ifft(sqrt(G).*fft(Einy));

% scale noise components by frequency-dependent gain S
Eampnx = ifft(sqrt(S).*fft(Nx));
Eampny = ifft(sqrt(S).*fft(Ny));

Eamp = [Eampsx + Eampnx; Eampsy + Eampny];

end