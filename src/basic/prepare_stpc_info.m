function STPCInfo = prepare_stpc_info(config)
% PREPARE_STPC_INFO  Build struct to be saved as mssa_Information.mat
%
% Usage:
%   STPCInfo = prepare_stpc_info(config)

    STPCInfo = struct();
    STPCInfo.Turning_number = config.turningNumber;
    STPCInfo.S_bou          = config.S_bound;
    STPCInfo.N_bou          = config.N_bound;
    STPCInfo.p_use          = config.p_use;

end