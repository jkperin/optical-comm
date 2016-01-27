% Output of the Equalizer
if Plots('Equalizer')
    scrsz = get(0,'ScreenSize');
    figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)])

    title ('Output of the equalizer')
    subplot(3,2,1); 
    plot(real(xhat(1,measind)),imag(xhat(1,measind)),'.'); 
    xlabel('\itx\rm-equ_{1i}'); ylabel('\itx\rm-equ_{1q}'); 
    axis image;

    subplot(3,2,2); 
    plot(real(xhat(2,measind)),imag(xhat(2,measind)),'.'); 
    xlabel('\itx\rm-equ_{2i}'); ylabel('\itx\rm-equ_{2q}'); 
    axis image;

    % equalizer adaptation
    subplot(3,2,[3,4]); 
    sqerr = abs(eps(1,:)).^2 + abs(eps(2,:)).^2;
    smoothmeas = 32;
    semilogy(1:Nsymb,AdEq.LMS.muVec,'-',AdEq.LMS.startdd*[1,1],[0.5*min(AdEq.LMS.muVec),2*max(AdEq.LMS.muVec)],'--');
    axis([0 Nsymb 0.5*min(AdEq.LMS.muVec) 2*max(AdEq.LMS.muVec)]);
    title('Equalizer step size parameter vs. time');
    xlabel('Symbol interval \itk'); ylabel('\mu(\itk\rm)');
    text(AdEq.LMS.startdd+4,0.5*max(AdEq.LMS.muVec),['Start decision-directed at \itk\rm = ' num2str(AdEq.LMS.startdd,4)])

    legend('Squared error |\epsilon|^2(\itk\rm)',...
        ['Moving average over ' num2str(smoothmeas,2) ' symbols'], ...
        'Location', 'SouthWest');

    subplot(3,2,[5,6]);            
    semilogy(1:Nsymb,sqerr,'g-',1:Nsymb,smooth(sqerr,smoothmeas),'b-',...
        [startmeas,Nsymb],MSE_meas(l,n)*[1,1],'k--',[1,Nsymb],MSE_AWGN_expected(l,n)*[1,1],'k:',...
        startmeas*[1,1],[0.5*min(sqerr),2*max(sqerr)],'k--');
    %axis([0 Nsymb 2*min(sqerr) 4*max(sqerr)]);
    axis([0 Nsymb 0.5*MSE_AWGN_expected(l,n) 4*max(sqerr)]);
    title(['Carrier ' num2str(m) ', ' num2str(oversamp) ' samples/symbol, equalizer squared error vs. time']);
    xlabel('Symbol interval \itk'); ylabel('|\epsilon|^2(\itk\rm)');
    legend('Squared error |\epsilon|^2(\itk\rm)',['Squared error, moving average over ' num2str(smoothmeas,2) ' symbols'],...
        ['Mean squared error, starting at \itk\rm = ' num2str(startmeas,4)],...
        'Mean squared error, approx. min. for ideal case');
    text(startmeas+20,0.5*max(sqerr),['Start measuring at \itk\rm = ' num2str(startmeas,4)]);
end
            
            
            