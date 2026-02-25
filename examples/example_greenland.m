% EXAMPLE_greenland  Example: Greenland case
%
% Usage:
%   Run in MATLAB:
%       examples/example_greenland
%
% Description:
%   This script builds a configuration struct "config" and calls
%   run_case(config) to complete the Slepian -> MSSA -> STPC workflow.
%
% Main Reference:
% Ref1: Harig, C., & Simons, F. J. (2012). Mapping GreenlandAs mass loss in space and time. Proceedings of the National Academy of Sciences, 109(49), 19934-19937.
% Ref2: Ma, Z., Fok, H. S., Tenzer, R., & Chen, J. (2024). A novel Slepian approach for determining mass-term sea level from GRACE over the South China Sea. International Journal of Applied Earth Observation and Geoinformation, 132, 104065.
% Ref3: Gauer, L. M., Chanard, K., & Fleitout, L. (2023). Data‐driven gap filling and spatio‐temporal filtering of the GRACE and GRACE‐FO records. Journal of Geophysical Research: Solid Earth, 128(5), e2022JB025561.
%
clear; clc;

tic

%% Build configuration (later you can load it from JSON, MAT file, etc.)
config = struct();

%% Set project root directory (assuming example_yangtze.m is under project_root)
config.ifilesRoot = 'D:\OneDrive\Slepian\Publication\release 3';
addpath(genpath(fullfile(config.ifilesRoot, 'src')));  % Add src to the MATLAB path

setenv('IFILES',config.ifilesRoot);  % define the route for Slepian processing

%% Code version
config.codeVersionBase = 'V2';

%% Region-related settings
% Make sure you have prepared the according code (e.g., 'yangtze.m') in
% /src/REGIONS_ADD/ (e.g., /src/REGIONS_ADD/yangtze.m)
% Notice that for oceanic applications (or specifically for leak-in
% simultion), another code for inner buffer zone should be defined. E.g.,
% SCSpTH.m (normal one) and SCSpTH_C.m (inner buffer one).
config.TH_ori        = 'greenland'; % the name of your study region 
config.areaName      = config.TH_ori; % the name of your study region
config.landOrOcean   = 'ice'; % it's for land or ocean? If for ocean, GAD component will be restored and inverted barameter (IB) correction will be calcualted (you can optionally add it to achieve mass sea level). 
config.c11cmn        = [-84.5 86.5 -4.5 55.5];  % regional boundary [lon_lower lat_upper lon_upper lat_lower]
% default resolution of results: 1 degree (this can not be changed)

%% GRACE/GRACE-FO Information
% Maximum spherical harmonic degree
config.Lwindow       = 60; % default degrees = 60
config.GIA           = 'Pelt17';

% Slepian data product
config.Institu_ver   = 'RL06';
% config.givInstitu    = ["CSR0324","JPL0324","GFZ0324","ITSG0324"]; % the name of your subfiles in GRACE/Originals/RL06(05)/*
config.givInstitu    = ["CSR0324","ITSG0324"]; % the name of your subfiles in GRACE/Originals/RL06(05)/*
% You should obey the rule of naming: [X..X]AABB.
% [X..X] is the abbreviation of your data (you name it).
% AA should represent the beginning year, and BB the ending year. 
% For example, CSR0323 means CSR GRACE and GRACE-FO data during 2003-2023.

%% STPC-related settings (significance level you want to get)
config.p_use         = [0.05, 0.1, 0.3]; % the significance level for consideration (e.g., 5%, 10%, 30%)

%% Code mode
% auto mode
config.buffer_deg    = 'Auto';   % how large is your buffer zone? select one (the number of your chosen one) in 'groupBuffer'.
config.M             = 'Auto'; 

% you can also define the groupBuffer if you like.
% NOTICE it is NEGATIVE for ocean, POSITIVE for land and ice.
config.groupBuffer   = [0,0.5,1,1.5];  % the range of groupBuffer your want to try

% customized mode
% config.buffer_deg    = 1;   % how large is your buffer zone? select one (the number of your chosen one) in 'groupBuffer'.
% % NOTICE it is NEGATIVE for ocean, POSITIVE for land and ice.
% config.groupBuffer   = [0,0.5,1,1.5]; % the range of groupBuffer your want to try
% config.M             = 120;

%% Do not change the following parameter if you don't know it.
%%
ensInstitu=[];
for i=1:numel(config.givInstitu)
    P_data=char(config.givInstitu(i));
    Pcenter=P_data(1:end-4);
    Ptime=P_data(end-3:end);
    ensInstitu=[ensInstitu Pcenter(1)];
end
ensInstitu=[ensInstitu Ptime];
config.ensInstitu    = ensInstitu; % the ensemble title according to the intitutions you used (that YOU NAME IT).
config.use_institu   = [config.givInstitu config.ensInstitu];

%% Slepian-related settings
% group of buffer zone ranges for usage. Positive number means outside the
% study region (normally for land), and negative number means inside the study region (normally for ocean).
config.S_choice      = 6; % the mode of selection of S for comparison (default mode 6 with fixed S = Max_S)
config.Max_S         = 50; % the maximum S only for plotting
config.Radius        = 500; % (km); if there is additional Gaussian filter
config.phi           = 0; % rotate coefficient
config.theta         = 0; % rotate coefficient
config.omega         = 0; % rotate coefficient
config.artificialMonths = 5*12+1:6*12; % artificial months for missing values interpolation and sensitivity analysis

%% MSSA-related settings for gap-filling
config.M_gap          = 13; % sliding window for gap filling. 13 months (Gauer et al., 2023)
config.N_gap          = 8; % cutoff number for gap filling. first 8 PCs

%% MSSA-related settings for reconstruction
config.turningNumber = 5; % the number of turning point for two lambda in STPC filter
config.S_bound       = 0.05; % the boundary restriction for lambda_S
config.N_bound       = 0; % the boundary restriction for lambda_N

%% Output directory (user should modify to their own local path)
config.resultDir = fullfile(config.ifilesRoot, ['Results_' config.codeVersionBase config.areaName]);
config.figureDir = fullfile(config.ifilesRoot, ['Figure_'  config.codeVersionBase config.areaName]);

% Options for figure plotting or not
config.Smooth       = ["None","Gaussian300km","DDK3","DDK5"];
config.redo         = false; % redo all the calculation
config.plotProcess  = true; % plot process figure
config.saveAddData  = true; % save complete data (otherwise only key data)
config.note         = ''; % any note for experiments

% 3. Run the case
results = run_case(config); 

disp([ config.areaName ' case finished.']);

toc
