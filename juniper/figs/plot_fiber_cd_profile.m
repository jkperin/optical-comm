clear, clc, close all

addpath ../f/
addpath ../data/cd/


filenames = {'CABLE-esmf_ocd.txt',...
    'CABLE-leaf_ocd.txt',...
    'CABLE-smf-3x25km_ocd.txt',...
    'CABLE-leaf-3x25km_ocd.txt',...
    'CABLE-smf-3km-1_ocd.txt',...
    'CABLE-smf-5km-1_ocd.txt',...
    'CABLE-smf-10km-1_ocd.txt',...
    'CABLE-smf-3km-2_ocd.txt',...
    'CABLE-smf-5km-2_ocd.txt',...
    'CABLE-smf-10km-2_ocd.txt'};

tags = {'dcf-smf', 'dcf-leaf', 'smf-3x25km', 'leaf-3x25km',...
    'smf-3km-1', 'smf-5km-1', 'smf-10km-1',...
    'smf-3km-2', 'smf-5km-2', 'smf-10km-2'};

CDprofiles = containers.Map();

for k = 1:length(filenames)
    [Wavelength, Delay, Dispersion, Slope] = read_cd_measurement_file(filenames{k});
    
    CDprofiles(tags{k}) = struct('wavelength', Wavelength, 'delay', Delay, 'D', Dispersion, 'S', Slope);
end

figure(1), hold on, box on
plot(CDprofiles('dcf-smf').wavelength, CDprofiles('smf-3x25km').D + CDprofiles('dcf-smf').D, 'Linewidth', 2)
plot(CDprofiles('dcf-leaf').wavelength, CDprofiles('leaf-3x25km').D + CDprofiles('dcf-leaf').D, 'Linewidth', 2)
plot(CDprofiles('dcf-smf').wavelength, CDprofiles('smf-3x25km').D + CDprofiles('dcf-smf').D + CDprofiles('smf-3km-1').D, 'Linewidth', 2)
plot(CDprofiles('dcf-smf').wavelength, CDprofiles('smf-3x25km').D + CDprofiles('dcf-smf').D + CDprofiles('smf-5km-1').D, 'Linewidth', 2)
plot(CDprofiles('dcf-smf').wavelength, CDprofiles('smf-3x25km').D + CDprofiles('dcf-smf').D + CDprofiles('smf-10km-1').D, 'Linewidth', 2)
plot(CDprofiles('dcf-smf').wavelength, CDprofiles('smf-3x25km').D + CDprofiles('dcf-smf').D + CDprofiles('smf-10km-1').D + CDprofiles('smf-3km-1').D, 'Linewidth', 2)
xlabel('Wavelength (nm)', 'FontSize', 12)
ylabel('Residual dispersion (ps/nm)', 'FontSize', 12)
legend('SMF and DCF', 'LEAF and DCF', 'SMF + DCF + 3km SMF', 'SMF + DCF + 5km SMF', 'SMF + DCF + 10km SMF', 'SMF + DCF + 13km SMF')

save('CDprofiles', 'CDprofiles')


figure, hold on, box on
plot(CDprofiles('smf-3x25km').wavelength, CDprofiles('smf-3x25km').D + CDprofiles('smf-10km-1').D + CDprofiles('smf-3km-1').D, 'b', 'Linewidth', 4)
xlabel('Wavelength (nm)', 'FontSize', 18)
ylabel('Dispersion (ps/nm)', 'FontSize', 18)
set(gca, 'FontSize', 18)
grid on
axis tight
axis([1530 1560 1380 1540])
set(gca, 'ytick', 1380:40:1540)

figure, hold on, box on
plot(CDprofiles('dcf-smf').wavelength, CDprofiles('dcf-smf').D, 'r', 'Linewidth', 4)
xlabel('Wavelength (nm)', 'FontSize', 18)
ylabel('Dispersion (ps/nm)', 'FontSize', 18)
set(gca, 'FontSize', 18)
grid on
axis tight
a = axis;
axis([1530 1560 -1540 -1380])
set(gca, 'ytick', -1540:40:-1380)


figure, hold on, box on
plot(CDprofiles('smf-3x25km').wavelength, CDprofiles('smf-3x25km').D + CDprofiles('dcf-smf').D + CDprofiles('smf-10km-1').D + CDprofiles('smf-3km-1').D, 'Color', [148,0,211]/255, 'Linewidth', 4)
xlabel('Wavelength (nm)', 'FontSize', 18)
ylabel('Dispersion (ps/nm)', 'FontSize', 18)
set(gca, 'FontSize', 18)
grid on
axis tight
axis([1530 1560 -20 20])
