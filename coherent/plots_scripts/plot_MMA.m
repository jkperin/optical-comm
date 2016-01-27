% Output of the Equalizer

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])

title ('Output of the equalizer')
subplot(2,2,1); 
plot(real(xhat(1,measind)),imag(xhat(1,measind)),'.'); 
xlabel('\itx\rm-equ_{1i}'); ylabel('\itx\rm-equ_{1q}'); 
axis image;
                        
subplot(2,2,2); 
plot(real(xhat(2,measind)),imag(xhat(2,measind)),'.'); 
xlabel('\itx\rm-equ_{2i}'); ylabel('\itx\rm-equ_{2q}'); 
axis image;
   

subplot(2,2,[3,4]); 
sqerr = abs(eps(1,:)).^2 + abs(eps(2,:)).^2;
smoothmeas = 32;
semilogy(1:Nsymb,sqerr,'g-', ...
    1:Nsymb,smooth(sqerr,smoothmeas),'b-', ...
    startmeas*[1,1],[0.5*min(sqerr),2*max(sqerr)],'k--');
title('Equalizer squared error vs. time')
axis([0 Nsymb .5*min(sqerr) 4*max(sqerr)])
xlabel('Symbol interval \itk'); ylabel('|\epsilon|^2(\itk\rm)');
legend('Squared error |\epsilon|^2(\itk\rm)',...
    ['Moving average over ' num2str(smoothmeas,2) ' symbols'], ...
    'Location', 'SouthWest');
text(startmeas+20,0.5*max(sqerr),['Start measuring at \itk\rm = ' num2str(startmeas,4)]);
