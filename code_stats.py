import os
from os import walk

local_path = "C:/Users/Joe/Dropbox/research/codes";

def file_len(fname):
    with open(fname, 'rb') as f:
        for i, l in enumerate(f):
            pass
    return i + 1


files = [];
for (dirpath, dirnames, filenames) in walk(local_path):
    for f in filenames:
        if ".git\\" in dirpath:
            continue;
        else:
            files.append(dirpath+"/"+f);


extensions = {};
for f in files:
    try:
    	p = f.index('.');
    except ValueError:
    	continue;

    try:
        extension = f[p:];
        lines = file_len(f);
        if extension in extensions:
                extensions[extension] += lines;
        else:
                extensions[extension] = lines;
                
    except Exception:
        print("Could not open file", f)

print("Extensions")
for e in extensions:
	print(e, extensions[e]);
