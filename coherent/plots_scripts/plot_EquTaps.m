figure()
if oversamp == 1.5 && strcmp(Equalization,'LMS')
    subplot(211)
    imagesc(abs(WTe));
    title(['Carrier ' num2set(m) ', equalizer coefficients for even \itk\rm |\itW\rm^{T}_{even}|'])
    colormap gray; axis image
    subplot(212)
    imagesc(abs(WTo));
    title(['Carrier ' num2set(m) ', equalizer coefficients for odd \itk\rm |\itW\rm^{T}_{odd}|'])
    colormap gray; axis image
else
    if strcmp(Equalization,'LMS')
        imagesc(abs(WTe));
    else 
        imagesc(abs(WT));
    end
    title(['Carrier ' num2str(m) ', equalizer coefficients |\itW\rm^{T}|'])
    colormap hot; axis image
end