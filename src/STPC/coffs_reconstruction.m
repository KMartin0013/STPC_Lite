function mssa_Sort_sn = coffs_reconstruction(MSSAcoffs_RCTest, MSSAcoffs_RC, fillM_deltacoffs, ...
                                           MSSAcoffs_rec, MSSAcoffs_evalues, S_ol, ...
                                           fill_nmonths, intit_num)
% processMSSAResults - Process MSSA analysis results, separate signal and noise components
%
% Input parameters:
%   MSSAcoffs_RCTest: Cell array of size S_ol×1, containing RC test results for each Slepian coefficient
%                     - Rows 1-2: Lilliefors and K-S test results (0/1)
%                     - Row 4: Frequency identification (1=annual, 2=semi-annual, <1=long-term, other=short-term)
%                     - Rows 5-7: Other parameters
%   MSSAcoffs_RC: Cell array of size S_ol×1, each element is an RC coefficient matrix
%                 - Dimensions: fill_nmonths × N_rec(ss) × intit_num
%   fillM_deltacoffs: 3D array, original Slepian coefficients
%                     - Dimensions: fill_nmonths × intit_num × S_ol
%   MSSAcoffs_rec: 3D array, reconstructed Slepian coefficients
%                  - Dimensions: fill_nmonths × intit_num × S_ol
%   MSSAcoffs_evalues: Eigenvalue matrix
%                      - Dimensions: N_rec(ss) × S_ol
%   S_ol: Scalar, number of Slepian coefficients
%   fill_nmonths: Scalar, length of time series (number of months)
%   intit_num: Scalar, number of initializations/iterations
%
% Output parameters:
%   mssa_Sort_sn: Structure array of size intit_num+1, containing the following fields:
%     - fillM_deltacoffs: Original Slepian coefficients
%     - fillM_reconcoffs: Reconstructed Slepian coefficients
%     - fillM_reconcoffs_RC: Cell array, coefficients for each RC
%     - fillM_reconcoffs_both: Sum of RCs that pass both tests (signal component)
%     - fillM_reconcoffs_noise: Sum of RCs that fail both tests (noise component)
%     - fillM_reconcoffs_resid: Residual component
%     - fillM_reconcoffs_eva: Eigenvalue statistics
%     - Frequency decomposition fields (annual, long-term, semi-annual, short-term, commented but can be enabled)
%
% Example:
%   mssa_Sort_sn = processMSSAResults(MSSAcoffs_RCTest, MSSAcoffs_RC, ...
%                                     fillM_deltacoffs, MSSAcoffs_rec, ...
%                                     MSSAcoffs_evalues, 10, 240, 5);

% Input parameter validation
if nargin < 8
    error('8 input parameters are required');
end

% Initialize output variable
mssa_Sort_sn = struct();

% Initialize intermediate variables
MSSAcoffs_RC_both = cell(S_ol, 5); % RCs that reject both tests, with frequency breakdown
MSSAcoffs_RC_noise = cell(S_ol, 5); % RCs that do not reject both tests, with frequency breakdown
MSSAcoffs_eva = []; % Eigenvalues summary statistics

RCFreDomi = cell(S_ol, 1);
RC_Seq = zeros(S_ol, 3);

% Initialize component matrices
MSSAcoffs_both = zeros(fill_nmonths, intit_num, S_ol); % Signal component
MSSAcoffs_noise = zeros(fill_nmonths, intit_num, S_ol); % Noise component
MSSAcoffs_resid = zeros(fill_nmonths, intit_num, S_ol); % Residual component

% Initialize frequency decomposition matrices (signal part)
MSSAcoffs_both_annual = zeros(fill_nmonths, intit_num, S_ol);
MSSAcoffs_both_long = zeros(fill_nmonths, intit_num, S_ol);
MSSAcoffs_both_semi = zeros(fill_nmonths, intit_num, S_ol);
MSSAcoffs_both_short = zeros(fill_nmonths, intit_num, S_ol);

% Initialize frequency decomposition matrices (noise part)
MSSAcoffs_noise_annual = zeros(fill_nmonths, intit_num, S_ol);
MSSAcoffs_noise_long = zeros(fill_nmonths, intit_num, S_ol);
MSSAcoffs_noise_semi = zeros(fill_nmonths, intit_num, S_ol);
MSSAcoffs_noise_short = zeros(fill_nmonths, intit_num, S_ol);

% Main processing loop: process each Slepian coefficient
MSSAcoffs_eva = zeros(S_ol, 11);
for ss = 1:S_ol
    % Get RC test results for current Slepian coefficient
    RCcoff_Test = MSSAcoffs_RCTest{ss};
    
    % Identify RCs that pass both tests (signal) and those that don't (noise)
    valid_both = find(RCcoff_Test(1, :) + RCcoff_Test(2, :) == 2);
    valid_noise = find(RCcoff_Test(1, :) + RCcoff_Test(2, :) < 2);
    
    % Record frequency information
    RCFreDomi{ss} = RCcoff_Test(4, :)';
    RC_Seq(ss, :) = mean(RCcoff_Test(5:7, :)');
    
    % Further classify signal RCs by frequency
    valid_both_annual = find(RCcoff_Test(1, :) + RCcoff_Test(2, :) == 2 & RCcoff_Test(4, :) == 1);
    valid_both_long = find(RCcoff_Test(1, :) + RCcoff_Test(2, :) == 2 & RCcoff_Test(4, :) < 1);
    valid_both_semi = find(RCcoff_Test(1, :) + RCcoff_Test(2, :) == 2 & RCcoff_Test(4, :) == 2);
    valid_both_short = setdiff(valid_both, [valid_both_annual, valid_both_long]); % Corrected version
    
    % Further classify noise RCs by frequency
    valid_noise_annual = find(RCcoff_Test(1, :) + RCcoff_Test(2, :) < 2 & RCcoff_Test(4, :) == 1);
    valid_noise_long = find(RCcoff_Test(1, :) + RCcoff_Test(2, :) < 2 & RCcoff_Test(4, :) < 1);
    valid_noise_semi = find(RCcoff_Test(1, :) + RCcoff_Test(2, :) < 2 & RCcoff_Test(4, :) == 2);
    valid_noise_short = setdiff(valid_noise, [valid_noise_annual, valid_noise_long]); % Corrected version
    
    % Store classification results
    MSSAcoffs_RC_both{ss, 1} = valid_both;
    MSSAcoffs_RC_both{ss, 2} = valid_both_long;
    MSSAcoffs_RC_both{ss, 3} = valid_both_annual;
    MSSAcoffs_RC_both{ss, 4} = valid_both_semi;
    MSSAcoffs_RC_both{ss, 5} = valid_both_short;
    
    MSSAcoffs_RC_noise{ss, 1} = valid_noise;
    MSSAcoffs_RC_noise{ss, 2} = valid_noise_long;
    MSSAcoffs_RC_noise{ss, 3} = valid_noise_annual;
    MSSAcoffs_RC_noise{ss, 4} = valid_noise_semi;
    MSSAcoffs_RC_noise{ss, 5} = valid_noise_short;
    
    % Get RC coefficients for current Slepian coefficient
    MSSAcoff_RC = MSSAcoffs_RC{ss};
    
    % Calculate signal component and its frequency decomposition
    MSSAcoffs_both(:, :, ss) = squeeze(sum(MSSAcoff_RC(:, valid_both, :), 2));
    MSSAcoffs_both_annual(:, :, ss) = squeeze(sum(MSSAcoff_RC(:, valid_both_annual, :), 2));
    MSSAcoffs_both_long(:, :, ss) = squeeze(sum(MSSAcoff_RC(:, valid_both_long, :), 2));
    MSSAcoffs_both_semi(:, :, ss) = squeeze(sum(MSSAcoff_RC(:, valid_both_semi, :), 2));
    MSSAcoffs_both_short(:, :, ss) = squeeze(sum(MSSAcoff_RC(:, valid_both_short, :), 2));
    
    % Calculate noise component and its frequency decomposition
    MSSAcoffs_noise(:, :, ss) = squeeze(sum(MSSAcoff_RC(:, valid_noise, :), 2));
    MSSAcoffs_noise_annual(:, :, ss) = squeeze(sum(MSSAcoff_RC(:, valid_noise_annual, :), 2));
    MSSAcoffs_noise_long(:, :, ss) = squeeze(sum(MSSAcoff_RC(:, valid_noise_long, :), 2));
    MSSAcoffs_noise_semi(:, :, ss) = squeeze(sum(MSSAcoff_RC(:, valid_noise_semi, :), 2));
    MSSAcoffs_noise_short(:, :, ss) = squeeze(sum(MSSAcoff_RC(:, valid_noise_short, :), 2));
    
    % Calculate residual component and recalculate noise component (for consistency)
    MSSAcoffs_resid(:, :, ss) = fillM_deltacoffs(:, :, ss) - squeeze(sum(MSSAcoff_RC, 2));
    MSSAcoffs_noise(:, :, ss) = squeeze(sum(MSSAcoff_RC, 2)) - squeeze(sum(MSSAcoff_RC(:, valid_both, :), 2));
    
    N_rec = size(MSSAcoff_RC, 2); % Number of RCs for current Slepian coefficient
    MSSAcoffs_eva(ss, :) = [...
        sum(MSSAcoffs_evalues(1:N_rec, ss), 1), ...                  % Sum of eigenvalues of all RCs
        sum(MSSAcoffs_evalues(valid_both, ss), 1), ...              % Sum of eigenvalues of signal RCs
        sum(MSSAcoffs_evalues(valid_both_long, ss), 1), ...         % Sum of eigenvalues of long-term signal RCs
        sum(MSSAcoffs_evalues(valid_both_annual, ss), 1), ...       % Sum of eigenvalues of annual signal RCs
        sum(MSSAcoffs_evalues(valid_both_semi, ss), 1), ...         % Sum of eigenvalues of semi-annual signal RCs
        sum(MSSAcoffs_evalues(valid_both_short, ss), 1), ...        % Sum of eigenvalues of short-term signal RCs
        sum(MSSAcoffs_evalues(valid_noise, ss), 1), ...             % Sum of eigenvalues of noise RCs
        sum(MSSAcoffs_evalues(valid_noise_long, ss), 1), ...        % Sum of eigenvalues of long-term noise RCs
        sum(MSSAcoffs_evalues(valid_noise_annual, ss), 1), ...      % Sum of eigenvalues of annual noise RCs
        sum(MSSAcoffs_evalues(valid_noise_semi, ss), 1), ...        % Sum of eigenvalues of semi-annual noise RCs
        sum(MSSAcoffs_evalues(valid_noise_short, ss), 1) ...        % Sum of eigenvalues of short-term noise RCs
    ];

end

% Create structure for each initialization result
for ins = 1:intit_num
    mssa_Sort_sn(ins).fillM_deltacoffs = squeeze(fillM_deltacoffs(:, ins, :));
    mssa_Sort_sn(ins).fillM_reconcoffs = squeeze(MSSAcoffs_rec(:, ins, :));
    
    % Store RC coefficients
    mssa_Sort_sn(ins).fillM_reconcoffs_RC = cell(S_ol, 1);
    for ss = 1:S_ol
        mssa_Sort_sn(ins).fillM_reconcoffs_RC{ss} = squeeze(MSSAcoffs_RC{ss}(:, :, ins));
    end
    
    % Store signal, noise and residual components
    mssa_Sort_sn(ins).fillM_reconcoffs_both = squeeze(MSSAcoffs_both(:, ins, :));
    mssa_Sort_sn(ins).fillM_reconcoffs_noise = squeeze(MSSAcoffs_noise(:, ins, :));
    mssa_Sort_sn(ins).fillM_reconcoffs_resid = squeeze(MSSAcoffs_resid(:, ins, :));
    
    % Store eigenvalue statistics
    mssa_Sort_sn(ins).fillM_reconcoffs_eva = MSSAcoffs_eva;
    
    % Optional: Enable frequency decomposition fields
    % mssa_Sort_sn(ins).fillM_reconcoffs_both_annual = squeeze(MSSAcoffs_both_annual(:, ins, :));
    % mssa_Sort_sn(ins).fillM_reconcoffs_both_long = squeeze(MSSAcoffs_both_long(:, ins, :));
    % mssa_Sort_sn(ins).fillM_reconcoffs_both_semi = squeeze(MSSAcoffs_both_semi(:, ins, :));
    % mssa_Sort_sn(ins).fillM_reconcoffs_both_short = squeeze(MSSAcoffs_both_short(:, ins, :));
end

% Create average results (average of all initializations)
mssa_Sort_sn(intit_num + 1).fillM_deltacoffs = squeeze(sum(fillM_deltacoffs, 2) / intit_num);
mssa_Sort_sn(intit_num + 1).fillM_reconcoffs = squeeze(sum(MSSAcoffs_rec, 2) / intit_num);

% Calculate average RC coefficients
mssa_Sort_sn(intit_num + 1).fillM_reconcoffs_RC = cell(S_ol, 1);
for ss = 1:S_ol
    RC_sum = zeros(size(mssa_Sort_sn(1).fillM_reconcoffs_RC{ss}));
    for ins = 1:intit_num
        RC_sum = RC_sum + mssa_Sort_sn(ins).fillM_reconcoffs_RC{ss};
    end
    mssa_Sort_sn(intit_num + 1).fillM_reconcoffs_RC{ss} = RC_sum / intit_num;
end

% Calculate average signal, noise and residual components
mssa_Sort_sn(intit_num + 1).fillM_reconcoffs_both = squeeze(sum(MSSAcoffs_both, 2) / intit_num);
mssa_Sort_sn(intit_num + 1).fillM_reconcoffs_noise = squeeze(sum(MSSAcoffs_noise, 2) / intit_num);
mssa_Sort_sn(intit_num + 1).fillM_reconcoffs_resid = squeeze(sum(MSSAcoffs_resid, 2) / intit_num);

% Store eigenvalue statistics and test results
mssa_Sort_sn(intit_num + 1).fillM_reconcoffs_eva = MSSAcoffs_eva;
mssa_Sort_sn(intit_num + 1).fillM_reconcoffs_RCTest = MSSAcoffs_RCTest;

end