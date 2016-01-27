% Output of the DD - phase Tracker
if (Plots('PhaseTracker'))
    scrsz = get(0,'ScreenSize');
    figure('OuterPosition',[1 1 scrsz(3) scrsz(4)/2])

    subplot(2,3,2); 
    plot(1:Nsymb,theta_pt(1,:))
    xlabel('Symbol interval \itk'); ylabel('Phase offset')
    axis([0 Nsymb min(theta_pt(1,:)-.5),max(theta_pt(1,:))+.5])

    subplot(2,3,5); 
    plot(1:Nsymb,theta_pt(2,:))
    xlabel('Symbol interval \itk'); ylabel('Phase offset')
    axis([0 Nsymb min(theta_pt(2,:))-.5,max(theta_pt(2,:))+.5])

    smoothmeas = 32;
    subplot(2,3,3);
    eps1_pt = abs(eps(1,:)).^2;
    semilogy(1:Nsymb,eps1_pt,'g-', ...
        1:Nsymb,smooth(eps1_pt,smoothmeas),'b-');
    title('Equalizer squared error vs. time')
    legend('Squared error |\epsilon|^2(\itk\rm)',...
        ['Moving average over ' num2str(smoothmeas,2) ' symbols'], ...
        'Location', 'SouthWest');
    axis([0 Nsymb .5*min(eps1_pt) 4*max(eps1_pt)])
    xlabel('Symbol interval \itk'); ylabel('|\epsilon|^2(\itk\rm)');

    subplot(2,3,6);
    eps2_pt = abs(eps(2,:)).^2;
    semilogy(1:Nsymb,eps2_pt,'g-', ...
        1:Nsymb,smooth(eps2_pt,smoothmeas),'b-');
    xlabel('Symbol interval \itk'); ylabel('|\epsilon|^2(\itk\rm)');

    subplot(2,3,1); 
    plot(real(xhat(1,measind)),imag(xhat(1,measind)),'.');
    xlabel('\itx\rm-equ_{1i}'); ylabel('\itx\rm-equ_{1q}'); 
    axis image;
    subplot(2,3,4); 
    plot(real(xhat(2,measind)),imag(xhat(2,measind)),'.');
    xlabel('\itx\rm-equ_{2i}'); ylabel('\itx\rm-equ_{2q}'); 
    axis image;
end