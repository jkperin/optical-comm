%function ParamOpt(inputfile)
%load(inputfile)
load('M1_Car6_OSNR10_EquCMA_nSC4.mat');   
% ---------- Parameter Range -----------
l_max=[40,60,80,100]
mu_CMA = 0.001:.002:.03;
mu_PT = 0.002:.003:.01;
LFilt = 4:7;
% --------- Performance Metric ---------
MSE_meas = zeros(length(mu_CMA),length(mu_PT),length(LFilt),length(l_max));
Qmeas = zeros(size(MSE_meas));
BER_from_Qmeas = zeros(size(MSE_meas)); 
BERmeas = zeros(size(MSE_meas)); 


Cpr.type        = 'DD';                                                    % carrier recovery type. 'dd' = decision-directed; 'nd' = non-data-aided
Cpr.Filter.Type = 'FIR';                                      

Cpr.Filter.f    = 0.95;
SNR=10;

% Decision-Directed Least Mean Squared Algorithm --------------------------
AdEq.LMS.LFilt    = 10; 
AdEq.LMS.muInit   = 0.02;                                          
AdEq.LMS.muFactor = .6;
AdEq.LMS.muInt    = 550;                                           
AdEq.LMS.muShifts = 10;
AdEq.LMS.muVec = AdEq.LMS.muInit*AdEq.LMS.muFactor.^(min((AdEq.LMS.muShifts+1),...      % vector of step size
    ceil((1:Nsymb)/AdEq.LMS.muInt))-1);     
% Training mode vs. decision-directed mode ----------------------------
AdEq.LMS.startdd = round(AdEq.LMS.muShifts/2)*AdEq.LMS.muInt;               % symbol index at which to stop using training symbols and change to decision-directed mode
AdEq.LMS.decwtvec = (1:Nsymb)<AdEq.LMS.startdd;                             % vector of decision weights (1 when use training symbols, 0 when use decision symbols)
startmeas = AdEq.LMS.muShifts*AdEq.LMS.muInt; 
Pt.startdd  = AdEq.LMS.startdd;    
Pt.decwtvec = (1:Nsymb)<Pt.startdd;                                         % vector of decision weights (1 when use training symbols, 0 when use decision symbols)

% -- Itr over different set of params --

for lmax = 1:length(l_max)
    for m1=1:length(mu_CMA) 
        for m2=1:length(mu_PT)
            for l=1:length(LFilt)
                clc
                display(sprintf('%.2f',m2/length(mu_PT)/length(mu_CMA)+(m1-1)/length(mu_CMA)));
                AdEq.CMA.mu    = mu_CMA(m1);
                AdEq.CMA.LFilt = LFilt(l);
                Pt.mu = mu_PT(m2);

                Cpr.Filter.Lmax = l_max(lmax);
                
                [xequ, eps_cma, WT] = CMA(tsamp, Ysep, Nsymb, oversamp, AdEq.CMA);
                [xcpr, ~, thetahat, Delta] = CPR(xequ,Cpr,PN,SNR);  
                [xD, xhat, eps, theta_pt] = PT(x,xcpr, Nsymb,Pt);
                %[xD, xhat, eps, theta_pt] = PT(x,xequ, Nsymb,Pt);

                measind = (startmeas:Nsymb)-Delta;                        % indices of symbols during which to measure                                                       
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
                MSE_meas(m1,m2,l,lmax) = mean(sqerr(measind));
                % Q factor (based on true error, starting at 'startmeas')
                Q1imeas = (mean(x1ihatmeas(x1imeas>0))-mean(x1ihatmeas(x1imeas<0)))/(std(eps1imeas(x1imeas<0))+std(eps1imeas(x1imeas>0)));
                Q1qmeas = (mean(x1qhatmeas(x1qmeas>0))-mean(x1qhatmeas(x1qmeas<0)))/(std(eps1qmeas(x1qmeas<0))+std(eps1qmeas(x1qmeas>0)));
                Q2imeas = (mean(x2ihatmeas(x2imeas>0))-mean(x2ihatmeas(x2imeas<0)))/(std(eps2imeas(x2imeas<0))+std(eps2imeas(x2imeas>0)));
                Q2qmeas = (mean(x2qhatmeas(x2qmeas>0))-mean(x2qhatmeas(x2qmeas<0)))/(std(eps2qmeas(x2qmeas<0))+std(eps2qmeas(x2qmeas>0)));            
                Qmeas(m1,m2,l,lmax) = (Q1imeas+Q1qmeas+Q2imeas+Q2qmeas)/4;
                % Bit-error ratio, average over all four tributaries, starting at 'startmeas', computed based on measured Q factors
                BER_from_Qmeas(m1,m2,l,lmax) = (0.5*erfc(Q1imeas/sqrt(2))+0.5*erfc(Q1qmeas/sqrt(2))+...
                    0.5*erfc(Q2imeas/sqrt(2))+0.5*erfc(Q2qmeas/sqrt(2)))/4;
                % measured
                num_bit_err_meas = sum(x1iDmeas~=x1imeas)+sum(x1qDmeas~=x1qmeas)+sum(x2iDmeas~=x2imeas)+sum(x2qDmeas~=x2qmeas);
                BERmeas(m1,m2,l,lmax) = num_bit_err_meas/4/numsymbmeas;

            end
        end
    end
end

for lmax = 1:length(l_max)
    figure
    for l = 1:length(LFilt)
        temp = Qmeas(:,:,l, lmax);
        ind = find(temp==max(max(temp)));
        r = mod(ind,length(mu_CMA));
        if(r==0)
            r=length(mu_CMA);
        end
        c = 1+(ind-r)/length(mu_CMA);        
        subplot(ceil(length(LFilt)/2),2,l)
        hold on
        imagesc(mu_PT,mu_CMA,temp)
        colormap autumn
        plot(mu_PT(c),mu_CMA(r),'--rs','LineWidth',2,...
             'MarkerEdgeColor','k',...
             'MarkerFaceColor','g',...
             'MarkerSize',10)
        dx = mean(diff(mu_PT));
        dy = mean(diff(mu_CMA));
        ratio=.99;
        A=mu_PT; B=mu_CMA; drawPatch;
        hold off
        xlabel('\mu_{PT}')
        ylabel('\mu_{CMA}')
        title(sprintf('# of taps: %d, %.3f',LFilt(l),temp(ind)));
        %axis([min(Factor) max(Factor) min(Init) max(Init)])
        axis tight
    end
end