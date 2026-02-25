function check_data(fullPathdata)

if exist(fullPathdata, 'dir') ~= 2
    error(['No such data: ' fullPathdata])
end

end