function save_waveform(x, filename)

if size(x, 1) < size(x, 2)
    x = x.';
end

try
    f = fopen(filename, 'w');
catch e
    fprintf('An error occured while trying to open file %s\n', filename)
    disp(e.message)
    return;
end

for k = 1:length(x)
    fprintf(f, '%d\n', x(k));
end

fclose(f);


