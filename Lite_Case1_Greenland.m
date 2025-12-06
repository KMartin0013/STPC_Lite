clc;clear;

%% select your workpath
setenv('IFILES','D:\OneDrive\Slepian\Publication\release 2');

addpath(fullfile(getenv('IFILES'),'slepian_abd-master_untouched'))
addpath(fullfile(getenv('IFILES'),'slepian_abd-master_untouched','REGIONS'))
addpath(fullfile(getenv('IFILES'),'slepian_abd-master_untouched','ICEREGIONS'))

addpath(fullfile(getenv('IFILES'),'STPC-Lite'))
%% Adjustable parameter (You can change it according to your needs.)

% data version
Code_Version='V1';

%%%% (case 1) for Greenland
% This is probably a case on land (assume that the change of outer ocean
% is negligible compared to inland.) Notice that oo=0 is a signature of
% land.

TH_ori='greenland';Area=TH_ori;
land_or_ocean='land';
ll=2;Code_Version=[Code_Version 'GL'];
c11cmn=[-84.5 86.5 -4.5 55.5]; % don't be so tough to yourself!

%%%%

%%%% (case 2) for South China Sea
% This is an example in the ocean.

% TH_ori='SCSpTH';Area=TH_ori;
% land_or_ocean='ocean';
% ll=[6];Code_Version=[Code_Version 'SCS'];
% c11cmn=[95.5 29.5 124.5 -4.5]; % don't be so tough to yourself!
%%%%

%%%% (case 3) for Yangtze
% This is a case on land (assume that the change of outer ocean
% is negligible compared to inland.)

% TH_ori='yangtze';Area=TH_ori;
% land_or_ocean='land';
% ll=3;Code_Version=[Code_Version 'Yang'];
% c11cmn=[85.5 39.5 124.5 20.5]; % don't be so tough to yourself!
%%%%

Lwindow=60; % the maximum level of SHC you want to use

ddir1=fullfile(getenv('IFILES'),['Results_' Code_Version]);
ddir2=fullfile(getenv('IFILES'),['Figure_' Code_Version]);

if ~exist(ddir1,'dir')
    mkdir(ddir1)
end
if ~exist(ddir2,'dir')
    mkdir(ddir2)
end

save('Basic_Information.mat','Code_Version','TH_ori','Area','Lwindow',...
    'land_or_ocean','c11cmn','ll','ddir1','ddir2')

%% Slepian part

% Decide on your choice of basis, depending on your region of interest
% and the bandwidth you want. Using GRACE2SLEPT, project the results of
% GRACE2PLMT into your chosen basis. We recommend that your basis is
% chosen based on a set of synthetic experiments which estimate the
% leakage/recovery tradeoffs.

% The name of product has a uniform format for ease of use. Please see more
% in 'grace2slept_m.m'.
Dataproduct={'JPL0323','RL06',Lwindow,land_or_ocean,'Pelt17'};

% this is a group of buffer zone by default
group_buffer=[0,0.5,1,1.5,-0.5,-1,-1.5,...
    -0.25,-0.75,-1.25,-2];
buffer_deg=group_buffer(ll); % degree of buffer zone

% Strategy of the choice of S. 
% Default is 0 for shannon number. 
% Case 1 for a larger number of N with most energy (i.e., 90%). 
% Case 2 is based on Case 1 but with additional Gaussian smooth when N is larger than shannon
% number. 
% Case 3 is all-smoothing version of case 1.
% Case 4 for most energy (i.e., 70%); Case 5 for most energy (i.e., 99%);
% Case 6 for S = 50; (make sure you have enough benchmark N for further calculation)
S_choice=6;Max_S=50; % A suitable SSF for Greenland, Yangtze River Basin and SCS

% Radius of additional Gaussian smooth
Radius=500; % Radius of additional Gaussian smooth

% this is by default for geographic areas
phi=0; theta=0; omega=0;

% Additional process for artificial months.
% If you're not sure you need it, odds are you don't need it. Please
% comment it out.
artificial_months=5*12+1:6*12;

if buffer_deg>=0
    buffer_str=num2str(buffer_deg);
    buffer_str(buffer_str=='.')='p';
else
    buffer_str=['neg' num2str(abs(buffer_deg))];
    buffer_str(buffer_str=='.')='p';
end

save('Slepian_Information','artificial_months',"omega","theta","phi",...
    "Radius","S_choice","Dataproduct","buffer_str","buffer_deg","group_buffer","Max_S")

% the Slepian procedure
Lite_Slepian_V1

disp('The time series of ''Code'' and ''Area-weighted'' should be generally consistent. Otherwise, check your ''c11cmn'' boundary of it contains the entire study region of interest (considering the buffer zones). Press any key to proceed.')
pause

%% Prepare all the Slepian coefficients for different data centers

disp('Before you do MSSA procedure, make sure you have run different data centers (if you tend to) for the Slepian part. Press any key to proceed.')
pause
close all

%% MSSA procedure

% How many institutions you want to do the M-SSA. Notice that the last
% component in 'use_institu' is the title of the ensemble output data than
% you can name it by yourself.
% example 1
% use_institu=["CSR0323","JPL0323","GFZ0323","ITSG0323","CJGI0323"];
% example 2
use_institu=["CSR0323","JPL0323","CJ0323"];
intit_num=numel(use_institu)-1;

% test for the lambda
% How many breaking points you want for the determination of S and N
% (for the lamda distribution)
Turning_number=5;

% the boundary restriction of S and N
% Normally we only use this to S, as the boundary constraints at the N
% exert a negligible influence on signals leaking externally.
% NOTICE: you had better not put the boundary too large.
S_bou=0.05;N_bou=0;

% the p values you want to use
SIG_use=[0.05,0.1,0.3];

save('MSSA_Information',"use_institu","intit_num","Turning_number","S_bou",...
    "N_bou","SIG_use")

Lite_MSSA_V1

%% STPC filter to determine the approxiamte N and S

Lite_STPC_Greenland_V1