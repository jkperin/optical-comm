%function ParamOpt(inputfile)
%load(inputfile)
load('M4_Car10_OSNR12_nSC1.mat');   
%load('M1_Car5_OSNR12.mat')
% ---------- Parameter Range -----------
Init = 0.01:.002:.09;
Factor = [.3:.1:.8];
Int  = [550];
Shift = [10];
AdEq.LFilt = 7;
% --------- Performance Metric ---------
MSE_meas = zeros(length(Init),length(Factor),length(Int),length(Shift));
Qmeas = zeros(size(MSE_meas));
BER_from_Qmeas = zeros(size(MSE_meas)); 
BERmeas = zeros(size(MSE_meas)); 

% -- Itr over different set of params --
for ini=1:length(Init) 
    for fa=1:length(Factor)
        for int=1:length(Int) 
            for sh = 1:length(Shift)
                clc
                display(sprintf('%.1f%%',100*((ini-1+(fa-1+((int-1)/length(Int)))/length(Factor))/length(Init))));
                AdEq.muInit   = Init(ini);                                                       
                AdEq.muFactor = Factor(fa);                         
                AdEq.muInt    = Int(int);
                AdEq.muShifts = Shift(sh);                                                          
                AdEq.muVec = AdEq.muInit*AdEq.muFactor.^(min((AdEq.muShifts+1),ceil((1:Nsymb)/AdEq.muInt))-1);      
                AdEq.startdd = round(AdEq.muShifts/2)*AdEq.muInt;                           
                AdEq.decwtvec = (1:Nsymb)<AdEq.startdd;                                     
                startmeas = AdEq.muShifts*AdEq.muInt;    
                [xhat, xD, eps, epsD, WTe, WTo] = LMS(x, tsamp, Ysep, Nsymb, oversamp, AdEq);
                % Performance metrics -----------------------------------------
                % Compute equalizer error, Q and BER
                measind = startmeas:Nsymb;                                      % indices of symbols during which to measure
                numsymbmeas = length(measind);                                  % number of symbols measured
                x1imeas = real(x(1,measind)); x1qmeas = imag(x(1,measind));             % correct symbols over indices measured
                x2imeas = real(x(2,measind)); x2qmeas = imag(x(2,measind));
                x1ihatmeas = real(xhat(1,measind)); x1qhatmeas = imag(xhat(1,measind)); % estimated symbols over indices measured
                x2ihatmeas = real(xhat(2,measind)); x2qhatmeas = imag(xhat(2,measind));
                x1iDmeas = real(xD(1,measind)); x1qDmeas = imag(xD(1,measind));         % decision symbols over indices measured
                x2iDmeas = real(xD(2,measind)); x2qDmeas = imag(xD(2,measind));
                eps1imeas = real(eps(1,measind)); eps1qmeas = imag(eps(1,measind));     % true error for indices measured
                eps2imeas = real(eps(2,measind)); eps2qmeas = imag(eps(2,measind));
                % True error and true mean-square error (starting at 'startmeas')
                sqerr = abs(eps(1,:)).^2 + abs(eps(2,:)).^2;
                MSE_meas(ini,fa,int,sh) = mean(sqerr(measind));
                % Q factor (based on true error, starting at 'startmeas')
                Q1imeas = (mean(x1ihatmeas(x1imeas>0))-mean(x1ihatmeas(x1imeas<0)))/(std(eps1imeas(x1imeas<0))+std(eps1imeas(x1imeas>0)));
                Q1qmeas = (mean(x1qhatmeas(x1qmeas>0))-mean(x1qhatmeas(x1qmeas<0)))/(std(eps1qmeas(x1qmeas<0))+std(eps1qmeas(x1qmeas>0)));
                Q2imeas = (mean(x2ihatmeas(x2imeas>0))-mean(x2ihatmeas(x2imeas<0)))/(std(eps2imeas(x2imeas<0))+std(eps2imeas(x2imeas>0)));
                Q2qmeas = (mean(x2qhatmeas(x2qmeas>0))-mean(x2qhatmeas(x2qmeas<0)))/(std(eps2qmeas(x2qmeas<0))+std(eps2qmeas(x2qmeas>0)));
                Qmeas(ini,fa,int,sh) = (Q1imeas+Q1qmeas+Q2imeas+Q2qmeas)/4;
                % Bit-error ratio, average over all four tributaries, starting at 'startmeas', computed based on measured Q factors
                BER_from_Qmeas(ini,fa,int,sh) = (0.5*erfc(Q1imeas/sqrt(2))+0.5*erfc(Q1qmeas/sqrt(2))+...
                0.5*erfc(Q2imeas/sqrt(2))+0.5*erfc(Q2qmeas/sqrt(2)))/4;
                % measured
                num_bit_err_meas = sum(x1iDmeas~=x1imeas)+sum(x1qDmeas~=x1qmeas)+sum(x2iDmeas~=x2imeas)+sum(x2qDmeas~=x2qmeas);
                BERmeas(ini,fa,int,sh) = num_bit_err_meas/4/numsymbmeas;
            end
        end
    end
end
figure()
for sh = 1:length(Shift)
    for int=1:length(Int)
        temp = Qmeas(:,:,int,sh);
        ind = find(temp==max(max(temp)));
        r= mod(ind,length(Init));
        if(r==0)
            r=length(Init);
        end
        c = 1+(ind-r)/length(Init);        
        subplot(length(Shift),length(Int),(sh-1)*length(Int)+int)
        hold on
        imagesc(Factor,Init,temp)
        colormap autumn
        plot(Factor(c),Init(r),'--rs','LineWidth',2,...
             'MarkerEdgeColor','k',...
             'MarkerFaceColor','g',...
             'MarkerSize',10)
        dx = mean(diff(Factor));
        dy = mean(diff(Init));
        ratio=.99;
        A=Factor; B=Init;
        drawPatch;
        hold off
        xlabel('\mu_{factor}')
        ylabel('\mu_{initial}')
        title(sprintf('Lfilt: %d, # of shifts: %d,Interval: %d, %.3f',AdEq.LFilt, Shift(sh),Int(int),temp(ind)));
        %axis([min(Factor) max(Factor) min(Init) max(Init)])
        axis tight
    end
end