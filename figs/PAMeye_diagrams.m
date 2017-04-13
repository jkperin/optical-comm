clear, clc, close all

addpath ../f/

% Initialize system parameters
Fs = 10000;
Rs = 100;
nSamps = Fs/Rs;
rolloff = 0.99;
M = 2;
hMod = comm.PAMModulator; % comm.QPSKModulator System object
hMod.ModulationOrder = M;

% Square root raised cosine filters. Apply gain to normalize passband
% filter gain to unity.
hTxFlt = comm.RaisedCosineTransmitFilter('RolloffFactor', rolloff, ...
    'OutputSamplesPerSymbol', nSamps, ...
    'FilterSpanInSymbols', 6, ...
    'Gain', 9.9121);
hRxFlt = comm.RaisedCosineReceiveFilter('RolloffFactor', rolloff, ...
    'InputSamplesPerSymbol', nSamps, ...
    'FilterSpanInSymbols', 6, ...
    'DecimationFactor', 1, ...
    'Gain', 0.1009);


% Generate modulated and pulse shaped signal
frameLen = 1000;
msgData = randi([0 M-1],frameLen,1);
msgSymbols = step(hMod, msgData);
msgTx = step(hTxFlt, msgSymbols);

% Create an eye diagram object
eyeObj = commscope.eyediagram(...
    'SamplingFrequency', Fs, ...
    'SamplesPerSymbol', nSamps, ...
    'OperationMode', 'Real Signal')

% Update the eye diagram object with the transmitted signal
update(eyeObj, 0.5*msgTx);

%% Channel
EsNo = 25; SNR = EsNo - 10*log10(nSamps);

% Create an comm.AWGNChannel System object and set its NoiseMethod property
% to 'Signal to noise ratio (SNR)'. Set the 'SignalPower' property to the
% calculated input signal power.
hChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)',...
  'SNR', SNR);
hChan.SignalPower = (msgTx' * msgTx)/ length(msgTx);
msgRx = step(hChan, msgTx);
msgRxMf = step(hRxFlt, msgRx);

% Discard the existing eye diagram data
reset(eyeObj);

% Update the eye diagram object with the received signal
update(eyeObj, msgRxMf);

eyeObj.PlotType = '2D Line';

eyeObj.NumberOfStoredTraces = 100;
reset(eyeObj);

frameLength = 200;
msgData = randi([0 M-1],frameLen,1);
msgSymbols = step(hMod, msgData);
msgTx = step(hTxFlt, msgSymbols);

hChan.SignalPower = (msgTx' * msgTx)/ length(msgTx);
msgRx = step(hChan, msgTx);
msgRxMf = step(hRxFlt, msgRx);

update(eyeObj, msgRxMf);
eyeObj.PlotTimeOffset = 0;

switch M
    case 2
        filename = 'ook-eyediagram.tex';
    case 4
        filename = '4PAM-eyediagram.tex';
end

matlab2tikz(filename)


