function check_softwares(fullPathdata)

if exist(fullPathdata, 'dir') ~= 7
    error(['No such data: ' fullPathdata])
end

end