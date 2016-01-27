%function ParamOpt(inputfile)
%load(inputfile)
load('M3_Car2_OSNR12_nSC1.mat');   
% ---------- Parameter Range -----------
PhaseNoise = 0;
mu_CMA = [6:1:16]*1e-3;
mu_LMS = [100:50:700]*1e-3;
LFilt = 7;
Factor = [.6];
Shifts = [14];

% --------- Performance Metric ---------
MSE_meas = zeros(length(mu_CMA),length(mu_LMS),length(Factor),length(Shifts));
Qmeas = zeros(size(MSE_meas));
BER_from_Qmeas = zeros(size(MSE_meas)); 
BERmeas = zeros(size(MSE_meas)); 


AdEq.CMA_DD.muInt      = 500;                                           

AdEq.CMA_DD.hardSwitch = 3000;

Pt.mu = .005;
% Training mode vs. decision-directed mode ----------------------------
Pt.startdd  = 2750;                                                           % symbol index at which to stop using training symbols and change to decision-directed mode
Pt.decwtvec = (1:Nsymb)<Pt.startdd;   

% -- Itr over different set of params --
for m1=1:length(mu_CMA)
    for m2=1:length(mu_LMS)
        for s=1:length(Factor)
            for l=1:length(Shifts)
                clc
                display(sprintf('%.2f',(m1+m2/length(mu_LMS))/length(mu_CMA)));
                
                AdEq.CMA_DD.LFilt      = LFilt;
                AdEq.CMA_DD.CMA_mu     = mu_CMA(m1);
                AdEq.CMA_DD.muInit     = mu_LMS(m2);
                AdEq.CMA_DD.muFactor   = Factor(s);                                        
                AdEq.CMA_DD.muShifts   = Shifts(l);  
                AdEq.CMA_DD.muVec   = [ AdEq.CMA_DD.CMA_mu*ones(1,AdEq.CMA_DD.hardSwitch),...
                    AdEq.CMA_DD.muInit*AdEq.CMA_DD.muFactor.^(min((AdEq.CMA_DD.muShifts+1),...      % vector of step size
                    ceil((AdEq.CMA_DD.hardSwitch+1:Nsymb)/AdEq.CMA_DD.muInt))-1)];  
                
                startmeas = max(Nsymb - 3000+1);
                
                [xhat, xD, eps, WT] = CMA_DD(tsamp, Ysep, Nsymb, oversamp, AdEq.CMA_DD);
                [xD, xhat, eps, theta_pt] = PT(x,xhat, Nsymb,Pt);
                measind = startmeas:Nsymb;
                
                % indices of symbols during which to measure                                                       
                % Performance metrics -----------------------------------------
                % Compute equalizer error, Q and BER
                numsymbmeas = length(measind);                                  % number of symbols measured
                x1imeas = real(x(1,measind)); x1qmeas = imag(x(1,measind));     % correct symbols over indices measured
                x2imeas = real(x(2,measind)); x2qmeas = imag(x(2,measind));
                % estimated symbols over indices measured
                x1ihatmeas = real(xhat(1,measind)); x1qhatmeas = imag(xhat(1,measind)); 
                x2ihatmeas = real(xhat(2,measind)); x2qhatmeas = imag(xhat(2,measind));          
                % decision symbols over indices measured                       
                x1iDmeas = real(xD(1,measind)); x1qDmeas = imag(xD(1,measind));       
                x2iDmeas = real(xD(2,measind)); x2qDmeas = imag(xD(2,measind));
                % true error for indices measured
                eps1imeas = real(eps(1,measind)); eps1qmeas = imag(eps(1,measind));     
                eps2imeas = real(eps(2,measind)); eps2qmeas = imag(eps(2,measind));
                sqerr = abs(eps(1,:)).^2 + abs(eps(2,:)).^2;
                MSE_meas(m1,m2,s,l) = mean(sqerr(measind));
                % Q factor (based on true error, starting at 'startmeas')
                Q1imeas = (mean(x1ihatmeas(x1imeas>0))-mean(x1ihatmeas(x1imeas<0)))/(std(eps1imeas(x1imeas<0))+std(eps1imeas(x1imeas>0)));
                Q1qmeas = (mean(x1qhatmeas(x1qmeas>0))-mean(x1qhatmeas(x1qmeas<0)))/(std(eps1qmeas(x1qmeas<0))+std(eps1qmeas(x1qmeas>0)));
                Q2imeas = (mean(x2ihatmeas(x2imeas>0))-mean(x2ihatmeas(x2imeas<0)))/(std(eps2imeas(x2imeas<0))+std(eps2imeas(x2imeas>0)));
                Q2qmeas = (mean(x2qhatmeas(x2qmeas>0))-mean(x2qhatmeas(x2qmeas<0)))/(std(eps2qmeas(x2qmeas<0))+std(eps2qmeas(x2qmeas>0)));            
                Qmeas(m1,m2,s,l) = (Q1imeas+Q1qmeas+Q2imeas+Q2qmeas)/4;
                % Bit-error ratio, average over all four tributaries, starting at 'startmeas', computed based on measured Q factors
                BER_from_Qmeas(m1,m2,s,l) = (0.5*erfc(Q1imeas/sqrt(2))+0.5*erfc(Q1qmeas/sqrt(2))+...
                    0.5*erfc(Q2imeas/sqrt(2))+0.5*erfc(Q2qmeas/sqrt(2)))/4;
                % measured
                num_bit_err_meas = sum(x1iDmeas~=x1imeas)+sum(x1qDmeas~=x1qmeas)+sum(x2iDmeas~=x2imeas)+sum(x2qDmeas~=x2qmeas);
                BERmeas(m1,m2,s,l) = num_bit_err_meas/4/numsymbmeas;

            end
        end
    end
end

for l=1:length(Shifts)
    figure
    for s = 1:length(Factor)
        temp = Qmeas(:,:,s,l);
        ind = find(temp==max(max(temp)));
        r = mod(ind,length(mu_CMA))
        if(r==0)
            r=length(mu_CMA);
        end
        c = 1+(ind-r)/length(mu_CMA)        
        subplot(ceil(length(Factor)/2),2,s)
        hold on
        imagesc(mu_LMS,mu_CMA,temp)
        colormap autumn
        plot(mu_LMS(c),mu_CMA(r),'--rs','LineWidth',2,...
             'MarkerEdgeColor','k',...
             'MarkerFaceColor','g',...
             'MarkerSize',10)
        dy = mean(diff(mu_CMA));
        dx = mean(diff(mu_LMS));
        ratio=.99;
        A=mu_LMS; B=mu_CMA; drawPatch;
        hold off
        xlabel('\mu_{LMS}')
        ylabel('\mu_{CMA}')
        title(sprintf('# of shifts: %d, Factor @ %.3f, %.3f',Shifts(l),Factor(s), temp(ind)));
        %axis([min(Factor) max(Factor) min(Init) max(Init)])
        axis tight
    end
end