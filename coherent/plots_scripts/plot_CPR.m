% Laser Phase Noise -------------------------------------------------------
if (Plots('CarrierPhaseNoise') && PhaseNoise == 1)
    scrsz = get(0,'ScreenSize');
    figure('OuterPosition',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])

    title ('Output of the equalizer')
    subplot(2,2,1); 
    plot(real(xcpr(1,measind)),imag(xcpr(1,measind)),'.'); 
    xlabel('\itx\rm-equ_{1i}'); ylabel('\itx\rm-equ_{1q}'); 
    axis image;

    subplot(2,2,2); 
    plot(real(xcpr(2,measind)),imag(xcpr(2,measind)),'.'); 
    xlabel('\itx\rm-equ_{2i}'); ylabel('\itx\rm-equ_{2q}'); 
    axis image;


    subplot(2,2,[3,4]); 
    offset  =  mean(Tx_PN) - mean(thetahat);
    plot((1:Nsymb*Nstep)/Nstep,Tx_PN,'-g',1:Nsymb,thetahat+offset,'-b');
    legend('Actual Carrier Phase \theta (\itt\rm)','Phase Estimate \theta-hat (\itt\rm)');
    xlabel('Symbol interval \itk'); ylabel('Carrier Phase \theta')
    axis([0 max(measind) min(Tx_PN)-1,max(Tx_PN)+1])
    title(['Carrier Phase Noise ( Phase offset ' num2str(offset,4) ' )'])
end

