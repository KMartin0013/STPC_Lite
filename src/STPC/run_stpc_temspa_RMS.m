function RMS_Results = run_stpc_temspa_RMS( ...
    Turning_number, turn_V_four, turn_MSSA_four, used_sta_SN, ...
    fill_months, in, c11cmn_area, STPC_p_SN_results)

% COMPUTE_MSSA_RMS
% This function performs all calculations related to MSSA/STPC results:
%   1. Compute averaged MSSA window size (turn_MSSA_four_ave)
%   2. Prepare grid (lon/lat mesh, area weights)
%   3. Check polygon mask (region of interest)
%   4. Compute temporal and spatial RMS for multiple "disp_data" cases
%
% INPUTS
%   Turning_number  : maximum S/N index (scalar)
%   turn_V_four     : function handle or matrix, V(S, station) -> window length
%   turn_MSSA_four  : cell array, each cell is an MSSA result matrix (nn x station)
%   used_sta_SN     : index of station (or station pair) being used
%   mssa_Sort       : struct array, mssa_Sort(1).fill_dates used for time axis
%   in        : 1D longitude / latitude vectors
%   XY              : Nx2 polygon vertices [lon, lat] for region
%   STPC_p_SN_results  
%                   : cell array {S,N}.EWH.* containing fields
%                     (fillM_MSSA, fillM_STPC, fillM_Grid_MSSA, etc.)
%
% OUTPUT
%   RMS_Results : struct with fields:
%       .S_range, .N_range
%       .turn_MSSA_four_ave
%       .lonlon, .latlat
%       .lonlon_matl, .latlat_matl
%       .c11cmn_area             : area weights (lat x lon)
%       .in                      : logical mask inside polygon
%       .t_ss_mean_rms           : [disp_data, time]
%       .t_nn_mean_rms           : [disp_data, time]
%       .t_ss_nn_rms             : [disp_data, S, N]
%       .s_ss_mean_rms           : [disp_data, N]
%       .s_nn_mean_rms           : [disp_data, S]
%       .s_ss_nn_rms             : [disp_data, S, N]
%       .t_ns_rms                : [disp_data, N]  (diagonal S=N)
%       .s_ns_rms                : [disp_data, N]  (diagonal S=N)
%       .fill_months             : time index (1:n_month)
%

%% Basic ranges
S_range = 1:Turning_number;
N_range = 1:Turning_number;

%% Compute the mean turn_MSSA_four_ave(S,N)
turn_MSSA_four_ave = zeros(Turning_number, Turning_number);

for ss = 1:Turning_number
    for nn = 1:Turning_number

        % Number of MSSA segments (V) for this S at given station
        V_ss = turn_V_four(ss, used_sta_SN);

        mm_sum = 0;
        for vv = 1:V_ss
            % turn_MSSA_four{vv} is assumed to be a matrix [N x station]
            mm_sum = mm_sum + turn_MSSA_four{vv}(nn, used_sta_SN);
        end

        % Average over all V segments
        turn_MSSA_four_ave(ss, nn) = mm_sum / V_ss;
    end
end

%% RMS calculations for different disp_data cases

use_ss = S_range;  % all S
use_nn = N_range;  % all N

% Pre-allocate containers (dim1 = disp_data = 1..5)
%   disp_data:
%       1: EWH_sn.fillM_MSSA
%       2: MSSA - STPC
%       3: STPC
%       4: EWH_sn.fillM_sig_STPC
%       5: EWH_sn.fillM_res_STPC

n_disp = 5;
n_time = numel(fill_months);
nS     = numel(S_range);
nN     = numel(N_range);

t_ss_mean_rms = zeros(n_disp, nN);
t_nn_mean_rms = zeros(n_disp, nS);
t_ss_nn_rms   = zeros(n_disp, nS, nN);

s_ss_mean_rms = zeros(n_disp, nN);
s_nn_mean_rms = zeros(n_disp, nS);
s_ss_nn_rms   = zeros(n_disp, nS, nN);

t_ns_rms = zeros(n_disp, nS);  % S=N
s_ns_rms = zeros(n_disp, nS);  % S=N

for disp_data = 1:n_disp

    % Matrix to accumulate means over S and N
    t_ss_mean = zeros(numel(use_nn), n_time); % index = N, time
    t_nn_mean = zeros(numel(use_ss), n_time); % index = S, time

    % For spatial RMS (area-weighted)
    s_ss_rms_mean = zeros(numel(use_nn), 1);  % average over S
    s_nn_rms_mean = zeros(numel(use_ss), 1);  % average over N

    % Loop over all S and N
    for ss = use_ss
        for nn = use_nn

            % Extract EWH structure for given (S,N)
            EWH_sn = STPC_p_SN_results{ss, nn}.EWH;

            % Select temporal and spatial fields based on disp_data
            switch disp_data
                case 1
                    t_res = EWH_sn.fillM_MSSA;
                    r_res = EWH_sn.fillM_Grid_MSSA;
                case 2
                    t_res = EWH_sn.fillM_MSSA - EWH_sn.fillM_STPC;
                    r_res = EWH_sn.fillM_Grid_MSSA - EWH_sn.fillM_Grid_STPC;
                case 3
                    t_res = EWH_sn.fillM_STPC;
                    r_res = EWH_sn.fillM_Grid_STPC;
                case 4
                    t_res = EWH_sn.fillM_sig_STPC;
                    r_res = EWH_sn.fillM_Grid_sig_STPC;
                case 5
                    t_res = EWH_sn.fillM_res_STPC;
                    r_res = EWH_sn.fillM_Grid_res_STPC;
            end

            % s_rms: spatial RMS at each grid cell (over time)
            % r_res is [time x lat x lon] or [lat x lon x time], adjust if needed.
            s_rms = squeeze(rms(r_res, 1));  % -> [lat x lon]

            % Temporal mean accumulation
            t_ss_mean(nn,:) = t_ss_mean(nn,:) + t_res;
            t_nn_mean(ss,:) = t_nn_mean(ss,:) + t_res;

            % Store full time series for S,N pair
            t_ss_nn_rms(disp_data, ss, nn) = rms(t_res);

            % Area-weighted spatial mean RMS inside polygon
            % s_rms(in) and area(in) are 1D arrays over masked region
            spatial_rms_area = sum(s_rms(in) .* c11cmn_area(in)) / sum(c11cmn_area(in));

            s_ss_rms_mean(nn,1) = s_ss_rms_mean(nn,1) + spatial_rms_area;
            s_nn_rms_mean(ss,1) = s_nn_rms_mean(ss,1) + spatial_rms_area;

            s_ss_nn_rms(disp_data, ss, nn) = spatial_rms_area;

        end
    end

    % Convert accumulated sums to RMS/mean values

    % Temporal RMS of average over S / N
    t_ss_mean_rms(disp_data,:) = rms( (t_ss_mean' / numel(use_ss)) );
    t_nn_mean_rms(disp_data,:) = rms( (t_nn_mean' / numel(use_nn)) );

    % Area-weighted spatial RMS averaged over S / N
    s_ss_mean_rms(disp_data,:) = s_ss_rms_mean / numel(use_ss);
    s_nn_mean_rms(disp_data,:) = s_nn_rms_mean / numel(use_nn);

    % Diagonal S=N contributions
    for ii = 1:numel(use_nn)
        s_ns_rms(disp_data,ii) = s_ss_nn_rms(disp_data, ii, ii);
        t_ns_rms(disp_data,ii) = t_ss_nn_rms(disp_data, ii, ii);
    end
end


%% Pack outputs into a struct

RMS_Results = struct();

RMS_Results.turn_MSSA_four_ave = turn_MSSA_four_ave;

RMS_Results.t_ss_mean_rms = t_ss_mean_rms;
RMS_Results.t_nn_mean_rms = t_nn_mean_rms;
RMS_Results.t_ss_nn_rms   = t_ss_nn_rms;

RMS_Results.s_ss_mean_rms = s_ss_mean_rms;
RMS_Results.s_nn_mean_rms = s_nn_mean_rms;
RMS_Results.s_ss_nn_rms   = s_ss_nn_rms;

RMS_Results.t_ns_rms = t_ns_rms;
RMS_Results.s_ns_rms = s_ns_rms;

end
