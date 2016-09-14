function newfilename = check_filename(filename)
%% Check whether file named 'filename' already exists. If yes, new file name will be filename_(counter).ext
% filename must contain extension

if exist(filename) == 2
    p1 = strfind(filename, '(');
    p2 = strfind(filename, ')');
    p3 = strfind(filename, '.');
    ext = filename(p3:end);
    if isempty(p1)
        newfilename = [filename(1:p3-1) '(1)' ext];
    else
        count = round(str2double(filename(p1+1:p2-1)))+1;
        newfilename = [filename(1:p1-1) '(' num2str(count) filename(p2:end)];
    end
    newfilename = check_filename(newfilename);
else
    newfilename = filename;
end