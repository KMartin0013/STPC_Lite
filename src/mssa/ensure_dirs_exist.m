function ensure_dirs_exist(Attach)
% ENSURE_DIRS_EXIST  Create all required directories if they do not exist.

if ~exist(Attach.txt_path_ALL, 'dir'); mkdir(Attach.txt_path_ALL); end
if ~exist(Attach.fig_path_ALL, 'dir'); mkdir(Attach.fig_path_ALL); end

for i = 1:numel(Attach.txt_path_each)
    p1 = char(Attach.txt_path_each(i));
    p2 = char(Attach.fig_path_each(i));
    if ~exist(p1, 'dir'); mkdir(p1); end
    if ~exist(p2, 'dir'); mkdir(p2); end
end

end