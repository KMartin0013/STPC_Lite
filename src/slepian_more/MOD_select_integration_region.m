function [TH, XY, BasinArea] = MOD_select_integration_region( ...
    land_or_ocean, TH_ori, TH_buffer, XY_ori, XY_buffer, ...
        BasinArea_buffer, BasinArea_ori)
% SELECT_INTEGRATION_REGION  Decide which polygon is used for integration.
%
% For land basins (except Greenland), use the original region (no buffer).
% For ocean basins, include buffer in the integration region.
% For regions with severe leakage-out effect (e.g., Greenland), also
% include buffer zones to maximize the restored signals

    % notice that this is TH_ori for land usage, while being
    % TH_buffer for ocean usage.
    if strcmp(land_or_ocean,'ocean')
        TH          = TH_buffer;
        XY          = XY_buffer;
        BasinArea   = BasinArea_buffer;

    else
        
        TH          = TH_ori;
        XY          = XY_ori;
        BasinArea   = BasinArea_ori;
    end

end