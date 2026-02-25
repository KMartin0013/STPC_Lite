function mssaInfo = prepare_mssa_info(config)
% PREPARE_MSSA_INFO  Build struct to be saved as mssa_Information.mat
%
% Usage:
%   basicInfo = prepare_basic_info(config)

    mssaInfo = struct();
    mssaInfo.use_institu    = config.use_institu;
    mssaInfo.intit_num      = numel(config.use_institu) - 1;
    mssaInfo.M_gap          = config.M_gap;
    mssaInfo.N_gap          = config.N_gap;
    mssaInfo.M              = config.M;

end