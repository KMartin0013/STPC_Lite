function create_paths(fullPath)

if exist(fullPath, 'dir') ~= 7
    mkdir(fullPath)
end

end