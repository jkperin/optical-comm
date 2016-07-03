function [Wavelength, Delay, Dispersion, Slope] = read_cd_measurement_file(filename)

file = fopen(filename, 'r');

s = fscanf(file, '%c');

dataStart = strfind(s, 'Wavelength');

if dataStart == -1
    disp('Invalid file');
end
% 

data = strsplit(s(dataStart:end), '\t');

Headers = data(1:5);

entry = 1;
N = (length(data)-1)/5;
Wavelength = zeros(N, 1);
Delay = zeros(N, 1);
Dispersion = zeros(N, 1);
D = zeros(N, 1);
Slope = zeros(N, 1);
for k = 6:5:length(data)-1
    Wavelength(entry) = str2double(data{k});
    Delay(entry) = str2double(data{k+1});
    Dispersion(entry) = str2double(data{k+2});
    D(entry) = str2double(data{k+3});
    Slope(entry) = str2double(data{k+4});
    entry = entry+1;
end

Wavelength = Wavelength(1:entry-1);
Delay = Delay(1:entry-1);
Dispersion = Dispersion(1:entry-1);
D = D(1:entry-1);
Slope = Slope(1:entry-1);

table(Wavelength, Delay, Dispersion, Slope)

fclose(file);