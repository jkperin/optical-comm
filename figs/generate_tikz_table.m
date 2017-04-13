function generate_tikz_table(x, y, filename)
%% Generate tikz table to be used in generating plots in latex

if isempty(x) || isempty(y) || length(x) ~= length(y)
    warning('File not created!\n x and y dimensions must be non-zero and must have the same length')
    return
end

try
    fileID = fopen(filename, 'w');
catch e
    warning('Error opening the file %f:\n%s\n', filename, e.message)
    return
end

fprintf(fileID, 'table[row sep=crcr]{\n');

for k = 1:length(x)
    fprintf(fileID,'\t%.8f %.8f \\\\\n', x(k), y(k));
end

fprintf(fileID, '};\n');
fclose(fileID);

