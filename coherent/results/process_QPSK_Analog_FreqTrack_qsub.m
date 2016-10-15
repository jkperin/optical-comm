%% Process data saved by QPSK_BER_qsub.m
clear, clc, close all

addpath ../
addpath ../f/
addpath ../DSP/
addpath ../../f/
addpath ../../apd/
addpath ../../soa/

% QPSK_Analog_FreqTrack_OPLL-logic_Npol=1_linewidth=200kHz_delay=250ps_freqstep=0

Rb = 2*112e9;
Rs = Rb/(4);
BERtarget = 1.8e-4;
CPR = {'costas', 'logic'};
delay = 250;
Modulator = 'SiPhotonics';
ModBW = 30;
linewidth = 200;
FreqStepMHz = 200;

R = 1;
q = 1.60217662e-19;

SNRdB2PrxdBm = @(SNRdB) 10*log10(10^(SNRdB/10)*2*q*Rs/(R*1e-3));
SNRdBref = SNRreq(BERtarget, 4, 'QAM');
PrefdBm = -35.0798; % SNRdB2PrxdBm(SNRdBref);

LineStyle = {'-', '--'};
Marker = {'o', 's', 'v'};
Color = {[51, 105, 232]/255, [153,153,155]/255, [255,127,0]/255};

tol = 1e-3;
for n = 1
    for Ncpr = 1
        for f = 1:length(FreqStepMHz)
            filename = sprintf('QPSK_Analog_FreqTrack_OPLL-%s_Npol=%d_linewidth=%dkHz_delay=%dps_freqstep=%d.mat',...
                CPR{n}, Ncpr, linewidth, delay, FreqStepMHz(f));
            try 
                S = load(filename, '-mat');
                
                % Wrapped difference between signal and VCO phase
                Ntaps = 101;
                phaseError = 2*pi*S.Rx.LO.freqOffset.*(1:S.sim.N)/S.sim.fs - S.Sf;
                phaseErrorf = filtfilt(ones(1, Ntaps)/ Ntaps, 1, phaseError); % smoothing
                wphaseErrorf = asin(sin(phaseError)); % wrap

                % Measure time until locking again
                dt = 1/S.sim.fs;
                fVCO = diff(filtfilt(ones(1, Ntaps)/ Ntaps, 1, S.Sf))/(2*pi*dt); % instantaneous frequency of VCO
                dphi = diff(phaseErrorf);
                adphi = abs(dphi);
                phiS = cumsum(adphi);
                phiS = phiS/phiS(end);
                plock = find(phiS >= 1-tol, 1, 'first');
                Nlock = plock - S.Tstep;

                % Count cycle slips
                [peaks, idx] = findpeaks(wphaseErrorf(S.Tstep:plock));
                Ncycle_slips = sum(peaks >= 0.99*pi/2);
                idx = idx(peaks >= 0.99*pi/2);
                peaks = peaks(peaks >= 0.9*pi/2);
    
                figure(1), clf, hold on
                plot(fVCO)   
                plot(S.Rx.LO.freqOffset, 'k')
                plot(idx + S.Tstep - 1, fVCO(idx + S.Tstep - 1), 'ob')
                a = axis;
                axis([S.Tstep-20 plock a(3:4)])
                pause(0.5)
                
                figure(2), clc, hold on
                plot(wphaseErrorf)
                plot(idx + S.Tstep - 1, wphaseErrorf(idx + S.Tstep - 1), 'ob')
                a = axis;
                axis([S.Tstep-20 plock a(3:4)])
                pause(0.5) 
%                 
%                 figure(1), clf
%                 subplot(211), box on, hold on
%                 plot(t*1e9, freqOffset/1e6, 'k')
%                 plot(t(1:end-1)*1e9, fVCO)
%                 plot(t(idx)*1e9, fVCO(idx), 'xk')
%                 a = axis;
%                 plot(1e9*t(plock)*[1 1],  a(3:4), ':k')
%                 xlabel('Time (ns)')
%                 ylabel('Frequency (MHz)')
%                 axis tight
%                 subplot(212), hold on, box on
%                 plot(t*1e9, dphi)
%                 plot(t(idx)*1e9, dphi(idx), 'xk')
%                 a = axis;
%                 plot(1e9*t(plock)*[1 1],  a(3:4), ':k')
%                 xlabel('Time (ns)')
%                 ylabel('Wrapped phase error (rad)')
%                 drawnow

                
                
                
            catch e
                throw(e)
                continue
            end
        end
    end
end

