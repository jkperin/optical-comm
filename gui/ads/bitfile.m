clear, clc, close all

fileName = 'bitfile';

bitStream = randi([0 1], [2^20 1]);

fid = fopen([fileName '.txt'], 'w');

for k = 1:length(bitStream)
    fprintf(fid, '%d\n', bitStream(k));
end

fclose(fid);

save(fileName, 'bitStream');
    