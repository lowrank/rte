subfolders = {'./class/tracer/private'};
compiled = 1;

for i = 1:size(subfolders, 2)
    f = subfolders{i};
    subf = dir(f);
    if sum([subf.isdir]) == 0
        disp(f)
        mkdir(f);
        compiled = 0;
    end
end

if compiled == 0
!make all
end

clear;

