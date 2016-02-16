function [Vout, symbolsTX] = QAM_SC_Tx(dataTX, Tx, sim)
%% Generates QAM signal to drive electro optical modulator
% QAM signal generated has levels pre-compensated and is filtered by
% transmitter (DAC) filter Tx.filt
% Inputs: 
% - DataTX = symbol stream [0, ..., M-1], where M is QAM order
% - Tx, must contain constellation order Tx.const.M and transmitter filter Tx.filt
% - sim = simulation parameters (f, fs, Mct)

% Constellation
M = sim.M; % QAM order
modObj = modem.qammod('M', M, 'SymbolOrder', 'Gray');

X = modulate(modObj, dataTX(1, :)); % pol X
Y = modulate(modObj, dataTX(2, :)); % pol Y

% 
symbolsTX = [X; Y]; % transmitted symbols. Used in data-directed detection

% extract I and Q in x and y
Xi = real(X);
Xq = imag(X);
Yi = real(Y);
Yq = imag(Y);

% Compensates nonlinearity of electro-optical modulator
Xi_pd = asin(Xi/(log2(M)-1)); % normalize so that Xi in [-1, 1]
Xq_pd = asin(Xq/(log2(M)-1));
Yi_pd = asin(Yi/(log2(M)-1));
Yq_pd = asin(Yq/(log2(M)-1));

% create rectangular waveforms for modulators
V1x = Tx.Mod.Vpi*reshape(repmat(Xi_pd,sim.Mct,1), 1, sim.N);
V2x = Tx.Mod.Vpi*reshape(repmat(Xq_pd,sim.Mct,1), 1, sim.N);
V1y = Tx.Mod.Vpi*reshape(repmat(Yi_pd,sim.Mct,1), 1, sim.N);
V2y = Tx.Mod.Vpi*reshape(repmat(Yq_pd,sim.Mct,1), 1, sim.N);

% filter drive waveforms for modulators txfilt.
% group delay of Tx.filt.H has already been removed
Htx = ifftshift(Tx.filt.H(sim.f/sim.fs)); % transmitter filter
Vix = real(ifft(fft(V1x).*Htx));
Vqx = real(ifft(fft(V2x).*Htx)); 
Viy = real(ifft(fft(V1y).*Htx)); 
Vqy = real(ifft(fft(V2y).*Htx));

% Build output
Vout = [Vix + 1j*Vqx; Viy + 1j*Vqy];
