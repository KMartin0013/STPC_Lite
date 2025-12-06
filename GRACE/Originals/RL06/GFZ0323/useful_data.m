used_data_dir=dir('./');

use_data_name= struct('name', {used_data_dir(1:length(used_data_dir)).name});

save('used_files','use_data_name')



