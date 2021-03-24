function [folders] = find_folders(inputPath)

%find all subdirectories
filesPath = dir(inputPath);
dirFlags  = [filesPath.isdir];
folders = filesPath(dirFlags);
folders(ismember( {folders.name}, {'.', '..'})) = []; %%remove gar

end