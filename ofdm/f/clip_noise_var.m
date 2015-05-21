%% Calculate normalized clipping noise variance at each subcarrier

function varDk = clip_noise_var(dncp, ofdm, sig2)

% Clipping noise
Nc = ofdm.Nc;
Npre_os = ofdm.Npre_os;

Nsymb = length(dncp)/(Nc + Npre_os);

dncp = reshape(dncp, Nc + Npre_os, Nsymb);              % reshape into matrix form
dn = dncp(Npre_os+1:end, :);                           % remove cyclic prefix (does not consider ISI)
Dk = fft(dn, Nc, 1)/Nc; 
 
varDk = var(Dk, 1, 2);
varDk = varDk(2:Nc/2)/sig2;
varDk = varDk.';

figure
stem(varDk, 'fill')
xlabel('Subcarrier')
ylabel('Normalized clipping noise variance')
 
