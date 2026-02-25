function [mssa_Sort_EWH_sn,c11cmn_area,in,on] = grids_reconstruction_EWHorMSL( ...
    mssa_Sort_coff_sn,intit_num,fill_nmonths,S_rec,S_ol,r_record, ...
    lat,lon,XY,BasinArea)

% -------------------------------------------------------------------------
% reconstruct_gridded_EWH_MSL
%
% INPUT:
%   mssa_Sort_coff_sn   : struct array containing MSSA-related coefficients
%                       (fields must include: fill_SHCfilldelta, fillM_deltacoffs,
%                       fillM_reconcoffs, fillM_reconcoffs_both,
%                       fillM_IBdeltacoffs, etc.)
%   intit_num      : number of MSSA iterations (scalar), loop will be 1:intit_num+1
%   fill_nmonths   : number of months in the time series (scalar)
%   S_rec          : number of Slepian functions used for reconstruction (scalar)
%   r_record       : precomputed Slepian basis on grid, size [S_rec, nlat, nlon]
%   lat, lon       : 1-D latitude and longitude vectors defining the grid
%   lonlon, latlat : 2-D longitude/latitude meshgrid for the study area
%   XY             : polygon vertices (in lon-lat) for the basin mask
%   BasinArea      : basin area in m^2
%   S_ol           : index of IB coefficient in the slepian expansion
%
% OUTPUT:
%   mssa_Sort_sn   : updated struct with gridded EWH/MSL and regional
%                    (basin-averaged) variables written into the fields:
%                    fillM_Grid_EWH, fillM_Grid_EWH_MSSA, ...
%                    fillM_Grid_MASS, fillM_Grid_MASS_MSSA, ...
%                    fillM_EWH, fillM_EWH_MSSA, ..., fillM_MASS, ...
%                    fillM_MSL, fillM_MSL_MSSA, ..., fillM_Grid_MSL, etc.
%   c11cmn_area    : grid-cell area matrix (same size as lon/lat grid)
%   in, on         : logical masks from check_polygon_in, indicating grid
%                    cells inside / on the polygon XY
%
% NOTE:
%   - All computation steps and parameter names follow the original script.
%   - Only encapsulation into a function and indentation were changed.
% -------------------------------------------------------------------------

%% Reconstruct the gridded EWH (MSL)
% Some criteria for the calculation is that:
% we have the overall Mass and the gridded EWH(OBP) through the slepian
% calculation, then we need to calculate:
% (1) basin-averaged EWH (or basin-averaged MSL): from overall Mass
% (2) gridded Mass (or gridded MSL): from gridded EWH
%
% In addition, there is another process (M-SSA). Then, we will first
% calculate the corresponding overall MASS and the gridded EWH(OBP), then
% other variables as mentioned above.
%
% Besides, for the original SHC product, only the overall Mass and the
% gridded EWH(OBP) will be provided.

fitwhat = [3 365.25 365.25/2]; % linear trend and an annual and a semiannual harmonic term

% [~,N_SHC] = size(mssa_Sort_sn(1).fill_SHCfilldelta); %#ok<NASGU>

[lonlon,latlat] = meshgrid(lon, lat);
[in,on] = check_polygon_in(XY,lonlon,latlat); %#ok<ASGLU>

% area-weighted should be used to more precisely calculate the regional
% basin-averaged variables
c11cmn_area = zeros(length(lat),length(lon));
for i = 1:length(lat)
    for j = 1:length(lon)
        c11cmn_area(i,j) = areaquad(lat(i)-0.5,lon(j)-0.5, ...
                                    lat(i)+0.5,lon(j)+0.5);
    end
end

%% M-SSA (Slepian) part
% we first calculate the overall Mass and the gridded EWH(OBP) through the
% slepian calculation, then we need to calculate:
% (1) basin-averaged EWH (or basin-averaged MSL): from overall Mass
% (2) gridded Mass (or gridded MSL): from gridded EWH
%
% do selection criteria for both K-S and lillietest tests to rejected noise
% this is for the basin-averaged Mass or EWH (MSL)

mssa_Sort_EWH_sn=struct();
for ins = 1:intit_num+1

    % the original coefficients without MSSA has already been calculated
    % for the spatial and temporal results.
    fillM_deltacoffs_ins      = mssa_Sort_coff_sn(ins).fillM_deltacoffs;
%   fill_signaldelta_ins      = mssa_Sort_sn(ins).fill_signaldelta;

    % (spatial)
    % the gap-filling part with MSSA
    fillM_reconcoffs_ins      = mssa_Sort_coff_sn(ins).fillM_reconcoffs;
    fillM_reconcoffs_both_ins = mssa_Sort_coff_sn(ins).fillM_reconcoffs_both;

    for i = 1:fill_nmonths
        sp_ewh_mssa              = 0;
        sp_ewh_mssa_reconst      = 0;
        sp_ewh_mssa_both_reconst = 0;
        sp_ewh_mssa_resid_reconst     = 0;
        sp_ewh_mssa_bothfail_reconst  = 0;

        for j = 1:S_rec
            % r is the j-th Slepian function on the grid
            r = squeeze(r_record(j,:,:));

            sp_ewh_mssa = sp_ewh_mssa + ...
                r * fillM_deltacoffs_ins(i,j)' / 1000 * 100; % kg/m^2 -> cm
            sp_ewh_mssa_reconst = sp_ewh_mssa_reconst + ...
                r * fillM_reconcoffs_ins(i,j)' / 1000 * 100; % kg/m^2 -> cm
            sp_ewh_mssa_both_reconst = sp_ewh_mssa_both_reconst + ...
                r * fillM_reconcoffs_both_ins(i,j)' / 1000 * 100; % kg/m^2 -> cm
            sp_ewh_mssa_resid_reconst = sp_ewh_mssa_resid_reconst + ...
                r * (fillM_deltacoffs_ins(i,j)' - fillM_reconcoffs_ins(i,j)') / 1000 * 100; % kg/m^2 -> cm
            sp_ewh_mssa_bothfail_reconst = sp_ewh_mssa_bothfail_reconst + ...
                r * (fillM_reconcoffs_ins(i,j)' - fillM_reconcoffs_both_ins(i,j)') / 1000 * 100; % kg/m^2 -> cm

        end

        % store gridded EWH (cm)
        mssa_Sort_EWH_sn(ins).fillM_Grid_EWH(i,:,:)               = sp_ewh_mssa;
        mssa_Sort_EWH_sn(ins).fillM_Grid_EWH_MSSA(i,:,:)          = sp_ewh_mssa_reconst;
        mssa_Sort_EWH_sn(ins).fillM_Grid_EWH_MSSA_both(i,:,:)     = sp_ewh_mssa_both_reconst;
        mssa_Sort_EWH_sn(ins).fillM_Grid_EWH_MSSA_resid(i,:,:)    = sp_ewh_mssa_resid_reconst;
        mssa_Sort_EWH_sn(ins).fillM_Grid_EWH_MSSA_bothfail(i,:,:) = sp_ewh_mssa_bothfail_reconst;

        % convert gridded EWH (cm) to gridded MASS (Gt)
        % factor: (EWH_cm/100) * area(m^2) * 4*pi*R^2 * density(kg/m^3) / (1000 kg/t) / (10^9 t/Gt)
        mssa_Sort_EWH_sn(ins).fillM_Grid_MASS(i,:,:)               = sp_ewh_mssa              /100 .* c11cmn_area * 4*pi*6371000^2 * 10^3 / 10^3 / 10^9;
        mssa_Sort_EWH_sn(ins).fillM_Grid_MASS_MSSA(i,:,:)          = sp_ewh_mssa_reconst      /100 .* c11cmn_area * 4*pi*6371000^2 * 10^3 / 10^3 / 10^9;
        mssa_Sort_EWH_sn(ins).fillM_Grid_MASS_MSSA_both(i,:,:)     = sp_ewh_mssa_both_reconst /100 .* c11cmn_area * 4*pi*6371000^2 * 10^3 / 10^3 / 10^9;
        mssa_Sort_EWH_sn(ins).fillM_Grid_MASS_MSSA_resid(i,:,:)    = sp_ewh_mssa_resid_reconst/100 .* c11cmn_area * 4*pi*6371000^2 * 10^3 / 10^3 / 10^9;
        mssa_Sort_EWH_sn(ins).fillM_Grid_MASS_MSSA_bothfail(i,:,:) = sp_ewh_mssa_bothfail_reconst/100 .* c11cmn_area * 4*pi*6371000^2 * 10^3 / 10^3 / 10^9;

        % regional (basin-averaged) EWH (cm)
        mssa_Sort_EWH_sn(ins).fillM_EWH(i)               = sum(sp_ewh_mssa(in)              .* c11cmn_area(in)) / sum(c11cmn_area(in));
        mssa_Sort_EWH_sn(ins).fillM_EWH_MSSA(i)          = sum(sp_ewh_mssa_reconst(in)      .* c11cmn_area(in)) / sum(c11cmn_area(in));
        mssa_Sort_EWH_sn(ins).fillM_EWH_MSSA_both(i)     = sum(sp_ewh_mssa_both_reconst(in) .* c11cmn_area(in)) / sum(c11cmn_area(in));
        mssa_Sort_EWH_sn(ins).fillM_EWH_MSSA_resid(i)    = sum(sp_ewh_mssa_resid_reconst(in).* c11cmn_area(in)) / sum(c11cmn_area(in));
        mssa_Sort_EWH_sn(ins).fillM_EWH_MSSA_bothfail(i) = sum(sp_ewh_mssa_bothfail_reconst(in).* c11cmn_area(in)) / sum(c11cmn_area(in));

        % regional MASS (Gt), from basin-averaged EWH (cm)
        mssa_Sort_EWH_sn(ins).fillM_MASS(i)               = mssa_Sort_EWH_sn(ins).fillM_EWH(i)               /100 * BasinArea * 10^3 / 10^3 / 10^9;
        mssa_Sort_EWH_sn(ins).fillM_MASS_MSSA(i)          = mssa_Sort_EWH_sn(ins).fillM_EWH_MSSA(i)          /100 * BasinArea * 10^3 / 10^3 / 10^9;
        mssa_Sort_EWH_sn(ins).fillM_MASS_MSSA_both(i)     = mssa_Sort_EWH_sn(ins).fillM_EWH_MSSA_both(i)     /100 * BasinArea * 10^3 / 10^3 / 10^9;
        mssa_Sort_EWH_sn(ins).fillM_MASS_MSSA_resid(i)    = mssa_Sort_EWH_sn(ins).fillM_EWH_MSSA_resid(i)    /100 * BasinArea * 10^3 / 10^3 / 10^9;
        mssa_Sort_EWH_sn(ins).fillM_MASS_MSSA_bothfail(i) = mssa_Sort_EWH_sn(ins).fillM_EWH_MSSA_bothfail(i) /100 * BasinArea * 10^3 / 10^3 / 10^9;
    end

    if S_ol > S_rec
        for i = 1:fill_nmonths

            % basin-averaged MSL (cm), from basin-averaged EWH (cm) and IB term
            mssa_Sort_EWH_sn(ins).fillM_MSL(i) = ...
                mssa_Sort_EWH_sn(ins).fillM_EWH(i) * 1000/1028 - ...   % sea water
                fillM_deltacoffs_ins(i,S_ol) * 1000/1028/10;       % IB: mmH2O (kg/m^2) -> cm

            mssa_Sort_EWH_sn(ins).fillM_MSL_MSSA(i) = ...
                mssa_Sort_EWH_sn(ins).fillM_EWH_MSSA(i) * 1000/1028 - ...
                fillM_reconcoffs_ins(i,S_ol) * 1000/1028/10;

            mssa_Sort_EWH_sn(ins).fillM_MSL_MSSA_both(i) = ...
                mssa_Sort_EWH_sn(ins).fillM_EWH_MSSA_both(i) * 1000/1028 - ...
                fillM_reconcoffs_both_ins(i,S_ol) * 1000/1028/10;

            mssa_Sort_EWH_sn(ins).fillM_MSL_MSSA_resid(i) = ...
                mssa_Sort_EWH_sn(ins).fillM_EWH_MSSA_resid(i) * 1000/1028 - ...
                (fillM_deltacoffs_ins(i,S_ol) - fillM_reconcoffs_ins(i,S_ol)) * 1000/1028/10;

            mssa_Sort_EWH_sn(ins).fillM_MSL_MSSA_bothfail(i) = ...
                mssa_Sort_EWH_sn(ins).fillM_EWH_MSSA_bothfail(i) * 1000/1028 - ...
                (fillM_reconcoffs_ins(i,S_ol) - fillM_reconcoffs_both_ins(i,S_ol)) * 1000/1028/10;

            % transfer EWH to OBP and then subtract IB to gain MSL (gridded)
            mssa_Sort_EWH_sn(ins).fillM_Grid_MSL(i,:,:) = ...
                mssa_Sort_EWH_sn(ins).fillM_Grid_EWH(i,:,:) * 1000/1028 - ...
                fillM_deltacoffs_ins(i,S_ol) * 1000/1028/10;

            mssa_Sort_EWH_sn(ins).fillM_Grid_MSL_MSSA(i,:,:) = ...
                mssa_Sort_EWH_sn(ins).fillM_Grid_EWH_MSSA(i,:,:) * 1000/1028 - ...
                fillM_reconcoffs_ins(i,S_ol) * 1000/1028/10;

            mssa_Sort_EWH_sn(ins).fillM_Grid_MSL_MSSA_both(i,:,:) = ...
                mssa_Sort_EWH_sn(ins).fillM_Grid_EWH_MSSA_both(i,:,:) * 1000/1028 - ...
                fillM_reconcoffs_both_ins(i,S_ol) * 1000/1028/10;

            mssa_Sort_EWH_sn(ins).fillM_Grid_MSL_MSSA_resid(i,:,:) = ...
                mssa_Sort_EWH_sn(ins).fillM_Grid_EWH_MSSA_resid(i,:,:) * 1000/1028 - ...
                (fillM_deltacoffs_ins(i,S_ol) - fillM_reconcoffs_ins(i,S_ol)) * 1000/1028/10;

            mssa_Sort_EWH_sn(ins).fillM_Grid_MSL_MSSA_bothfail(i,:,:) = ...
                mssa_Sort_EWH_sn(ins).fillM_Grid_EWH_MSSA_bothfail(i,:,:) * 1000/1028 - ...
                (fillM_reconcoffs_ins(i,S_ol) - fillM_reconcoffs_both_ins(i,S_ol)) * 1000/1028/10;

        end
    end
end

end
