function [XY_buffer, BasinArea_buffer, XY_ori, BasinArea_ori, Earth_radius] = ...
    compute_regions(TH_ori, buffer_deg)
% COMPUTE_REGIONS  Compute buffered and original polygons and their areas.
%
% Outputs:
%   XY_buffer        - [lon, lat] polygon with buffer
%   BasinArea_buffer - area [m^2] of buffered region
%   XY_ori           - original [lon, lat] polygon (no buffer)
%   BasinArea_ori    - area [m^2] of original region
%   Earth_radius     - radius used for area computation

    % Earth radius from fralmanac, consistent with GRACE processing
    Earth_radius = fralmanac('a_EGM96','Earth');

    % Buffered region
    XY_buffer = eval(sprintf('%s(%i,%f)', TH_ori, 0, buffer_deg));
    BasinArea_buffer = spharea(XY_buffer) * 4 * pi * Earth_radius^2;

    % Original region
    XY_ori = eval(sprintf('%s(%i,%f)', TH_ori, 0, 0));
    BasinArea_ori = spharea(XY_ori) * 4 * pi * Earth_radius^2;
end