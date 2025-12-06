% These series of Figures are used for JGR.
% Code version: 5.
% Code Note:
%   Figure 1:
%   Main features:
%       (1) show the global SHC grid with north-south striping
%       (2) show the study area
%       (3) show the slepian basis
%
% Tested on 9.13.0.2049777 (R2022b)
% Last modified by zhongtian.ma@connect.polyu.hk, 2/3/2024

% Main Reference:
% Ref1: Harig, C., & Simons, F. J. (2012). Mapping GreenlandAs mass loss in space and time. Proceedings of the National Academy of Sciences, 109(49), 19934-19937.
% Ref2: Ma, Z., Fok, H. S., Tenzer, R., & Chen, J. (2024). A novel Slepian approach for determining mass-term sea level from GRACE over the South China Sea. International Journal of Applied Earth Observation and Geoinformation, 132, 104065.
% Ref3: Gauer, L. M., Chanard, K., & Fleitout, L. (2023). Data‐driven gap filling and spatio‐temporal filtering of the GRACE and GRACE‐FO records. Journal of Geophysical Research: Solid Earth, 128(5), e2022JB025561.

clc;clear;

%% Adjustable parameter (You can change it according to your needs.)

load('Basic_Information.mat')
load('Slepian_Information.mat');
load('MSSA_Information.mat')

%% Main code
month_str=["January","February","March","April","May","June","July",...
    "August","September","October","November","December"];

threshold_S=[S_bou, 1-S_bou]; % restriction on extremely high or low lambda
threshold_N=[N_bou, 1-N_bou];

%%%% (case 1) for Greenland
% This is a case on land (assuming that the change of outer ocean
% is negligible compared to inland).

% parrent_path3=Code_Version;
% 
% vv=char(Dataproduct(2));
% suffix_area3=[vv(3:4) '_' Area];
% 
% suffix_XY_use3=[num2str(Lwindow) '_' buffer_str '_' num2str(Radius)];
%%%%

%%%% (case 2) for South China Sea
% This is an example in the ocean, assuming that the trend of the
% outer land is one-fifth that of the ocean, and the amplitude is
% five times that of the ocean.

parrent_path2=Code_Version;

vv=char(Dataproduct(2));
suffix_area2=[vv(3:4) '_' Area];

suffix_XY_use2=[num2str(Lwindow) '_' buffer_str '_' num2str(Radius)];
%%%%

%%%% (case 3) for Yangtze River basin
% This is a case on land (assuming that the change of outer ocean
% is negligible compared to inland).

% parrent_path1=[Code_Version 'Yang'];

% suffix_area1='06_yangtze';

% degree 1 buffer zone
% suffix_XY_use1='60_1_500';
%%%%

if S_bou==0 && N_bou==0
    suffix_bou=[];
elseif S_bou==0 
    suffix_bou=['_Hn' num2str(N_bou)];
elseif N_bou==0
    suffix_bou=['_Hs' num2str(S_bou)];
else
    suffix_bou=['_Hs' num2str(S_bou) 'n' num2str(N_bou) ];
end

% suffix_str1=[suffix_area1 '_' suffix_XY_use1  '_S' num2str(Max_S)];
% 
% % please rewrite the name of this path accordingly to the output of 'V5r1.m'
% Attach_ALL1=[char(use_institu(intit_num+1)) '_' suffix_str1 suffix_bou];
% txt_path_ALL1=fullfile(getenv('IFILES'),['Text_' parrent_path1],Attach_ALL1);
% fig_path_ALL1=fullfile(getenv('IFILES'),['Figure_' parrent_path1],Attach_ALL1);
% 
% XY_path1=fullfile(getenv('IFILES'),['Results_' parrent_path1],['XY_' ...
%     char(use_institu(1)) '_' suffix_area1 '_' suffix_XY_use1]);
% XY1=load(XY_path1);
% 
suffix_str2=[suffix_area2 '_' suffix_XY_use2  '_S' num2str(Max_S)];

% please rewrite the name of this path accordingly to the output of 'V5r1.m'
Attach_ALL2=[char(use_institu(intit_num+1)) '_' suffix_str2 suffix_bou];
txt_path_ALL2=fullfile(getenv('IFILES'),['Text_' parrent_path2],Attach_ALL2);
fig_path_ALL2=fullfile(getenv('IFILES'),['Figure_' parrent_path2],Attach_ALL2);

XY_path2=fullfile(getenv('IFILES'),['Results_' parrent_path2],['XY_' ...
    char(use_institu(1)) '_' suffix_area2 '_' suffix_XY_use2]);
XY2=load(XY_path2);

% suffix_str3=[suffix_area3 '_' suffix_XY_use3  '_S' num2str(Max_S)];
% 
% % please rewrite the name of this path accordingly to the output of 'V5r1.m'
% Attach_ALL3=[char(use_institu(intit_num+1)) '_' suffix_str3 suffix_bou];
% txt_path_ALL3=fullfile(getenv('IFILES'),['Text_' parrent_path3],Attach_ALL3);
% fig_path_ALL3=fullfile(getenv('IFILES'),['Figure_' parrent_path3],Attach_ALL3);
% 
% XY_path3=fullfile(getenv('IFILES'),['Results_' parrent_path3],['XY_' ...
%     char(use_institu(1)) '_' suffix_area3 '_' suffix_XY_use3 '.mat']);
% XY3=load(XY_path3);

%% basic information

%% for the first region
%%
for SIG=1:3
    % SIG=1;

    load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 'c11.mat']));

    % load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 '_TP_S7_N7_M96_compare90.mat']),...
    %     'turn_V_four','turn_MSSA_four','MaxNumChanges_MSSA','MaxNumChanges_V','XY','c11cmn_combine','c11cmn');
    load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 '_TP' num2str(Turning_number) ...
        '_S' num2str(Turning_number) '_N' num2str(Turning_number) '_M' num2str(M_rec) ...
        '_sig' num2str(SIG_use(SIG)) '_compare.mat']),...
        'turn_V_four','turn_MSSA_four','MaxNumChanges_MSSA','MaxNumChanges_V','XY',...
        'c11cmn_combine','c11cmn','M_rec','Noise_SigLev','used_sta_V','used_sta_MSSA');

    turn_V2_four=turn_V_four;turn_MSSA2_four=turn_MSSA_four;

    c11cmn2=c11cmn;
    S_range2=1:MaxNumChanges_V;N_range2=1:MaxNumChanges_MSSA;

    % Noise_SigLev=0.1;
    % used_sta_V=4;used_sta_MSSA=4;

    if MaxNumChanges_V==MaxNumChanges_MSSA
        MNC=MaxNumChanges_V;
    else
        error('For now, the number of turning points of V and MSSA should be the same.')
    end

    % calculate the mean M
    turn_MSSA2_four_ave=zeros(MNC,MNC);
    for ss=1:MNC
        for nn=1:MNC
            V_ss=turn_V2_four(ss,used_sta_V);

            mm_sum=0;
            for vv=1:V_ss
                mm_sum=mm_sum+turn_MSSA2_four{vv}(nn,used_sta_MSSA);
            end
            turn_MSSA2_four_ave(ss,nn)=mm_sum/V_ss;

        end
    end

    data_year_beg=A5_Centers_Sort(1).data_year_beg;
    fill_dates=A5_Centers_Sort(1).fill_dates;fill_months=1:numel(fill_dates);

    % V=A5_Centers_Sort(1).V;

    % c11cmn12=c11cmn_combine;
    % time wasting
    [r_example_combine12,lon_combine12,lat_combine12,Plm,degres]=plm2xyz(CC{1},1,c11cmn_combine,60);
    %     c11cmn2=[0.5 89.5 359.5 -89.5];
    %     [r_example2,lon,lat,Plm,degres]=plm2xyz(CC{1},1,c11cmn2,60);
    %         c11cmn3=[0.5 89.5 60.5 -49.5];
    %     [r_example3,lon,lat,Plm,degres]=plm2xyz(CC{1},1,c11cmn3,60);
    [lonlon_combine12,latlat_combine12]=meshgrid(lon_combine12,lat_combine12);

    lon_combine12_matl=lon_combine12-0.5;
    lat_combine12_matl=lat_combine12+0.5;
    [lonlon_combine12_matl,latlat_combine12_matl]=meshgrid(lon_combine12_matl,lat_combine12_matl);

    [r_example_c11cmn2,lon_c11cmn2,lat_c11cmn2,Plm,degres]=plm2xyz(CC{1},1,c11cmn2,60);
    %     c11cmn2=[0.5 89.5 359.5 -89.5];
    %     [r_example2,lon,lat,Plm,degres]=plm2xyz(CC{1},1,c11cmn2,60);
    %         c11cmn3=[0.5 89.5 60.5 -49.5];
    %     [r_example3,lon,lat,Plm,degres]=plm2xyz(CC{1},1,c11cmn3,60);
    [lonlon_c11cmn2,latlat_c11cmn2]=meshgrid(lon_c11cmn2,lat_c11cmn2);

    lon_c11cmn2_matl=lon_c11cmn2-0.5;
    lat_c11cmn2_matl=lat_c11cmn2+0.5;
    [lonlon_c11cmn2_matl,latlat_c11cmn2_matl]=meshgrid(lon_c11cmn2_matl,lat_c11cmn2_matl);

    c11cmn_comTo2_lat=(c11cmn_combine(2)-c11cmn2(2)+1):(c11cmn_combine(2)-c11cmn2(4)+1);
    c11cmn_comTo2_lon=(c11cmn2(1)-c11cmn_combine(1)+1):(c11cmn2(3)-c11cmn_combine(1)+1);

    [in2, on2]=check_polygon_in(XY,lonlon_combine12,latlat_combine12);

    % area-weighted should be used to more precisely calculate the regional
    % basin-averaged variables
    for i=1:length(lat_combine12)
        for j=1:length(lon_combine12)
            c11cmn12_area(i,j)=areaquad(lat_combine12(i)-0.5,lon_combine12(j)-0.5,lat_combine12(i)+0.5,lon_combine12(j)+0.5);
        end
    end

    ins=1;
    %% temporal STD (mean)

    use_ss=S_range2;
    use_nn=N_range2;

    t2_ss_mean_std=[];t2_ss_mean_rms=[];
    t2_nn_mean_std=[];t2_nn_mean_rms=[];
    t2_ss_nn_std=[];t2_ss_nn_rms=[];

    s2_ss_mean_std=[];s2_ss_mean_rms=[];
    s2_nn_mean_std=[];s2_nn_mean_rms=[];
    s2_ss_nn_std=[];s2_ss_nn_rms=[];
    for disp_data=1:5

        t_ss_mean=zeros(numel(use_nn),numel(fill_months));
        t_nn_mean=zeros(numel(use_ss),numel(fill_months));
        t_ss_nn=zeros(numel(use_ss),numel(use_nn),numel(fill_months));

        s_ss_std_mean=zeros(numel(use_nn),1);
        s_nn_std_mean=zeros(numel(use_ss),1);
        s_ss_rms_mean=zeros(numel(use_nn),1);
        s_nn_rms_mean=zeros(numel(use_ss),1);

        for ss=use_ss

            %         SSF=S_range2(ss);
            % choose nexttile or subplot (based on your matlab version containing tiledlayout or not)
            %             nexttile % nexttile method
            %         subplot(a,b,SSF)  % subplot method

            for nn=use_nn

                %             load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 '_TP_S' num2str(ss) '_N' num2str(nn) '_M60_compare.mat']),...
                %                 'A5_Centers_Compare','A5_Centers_STD_Compare');
                load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 '_TP' num2str(MNC) ...
                    '_S' num2str(ss) '_N' num2str(nn) '_M' num2str(M_rec) '_sig' num2str(Noise_SigLev) ...
                    '_compare.mat']),...
                    'A5_Centers_Compare','A5_Centers_STD_Compare','A5_Centers_RMS_Compare');

                switch disp_data
                    case 1
                        t_res=A5_Centers_Compare(ins).fill_EWH_nonoisecorrection;
                        s_std=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_uptoS_STD(end,:,:));
                        s_rms=squeeze(A5_Centers_RMS_Compare(ins).fill_Grid_EWH_MSSA_uptoS_RMS(end,:,:));
                        fig_name=['fill_EWH_nonoisecorrection'];
                    case 2
                        t_res=A5_Centers_Compare(ins).fill_EWH_noisecorrection;
                        s_std=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_bothfail_uptoS_STD(end,:,:));
                        s_rms=squeeze(A5_Centers_RMS_Compare(ins).fill_Grid_EWH_MSSA_bothfail_uptoS_RMS(end,:,:));
                        fig_name=['fill_EWH_noisecorrection'];
                    case 3
                        t_res=A5_Centers_Compare(ins).fill_EWH;
                        s_std=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_both_uptoS_STD(end,:,:));
                        s_rms=squeeze(A5_Centers_RMS_Compare(ins).fill_Grid_EWH_MSSA_both_uptoS_RMS(end,:,:));
                        fig_name=['fill_EWH'];
                    case 4
                        t_res=A5_Centers_Compare(ins).fill_EWHsig;
                        s_std=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWHsig_MSSA_both_uptoS_STD(end,:,:));
                        s_rms=squeeze(A5_Centers_RMS_Compare(ins).fill_Grid_EWHsig_MSSA_both_uptoS_RMS(end,:,:));
                        fig_name=['fill_EWHsig'];
                    case 5
                        t_res=A5_Centers_Compare(ins).fill_EWHres;
                        s_std=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWHres_MSSA_both_uptoS_STD(end,:,:));
                        s_rms=squeeze(A5_Centers_RMS_Compare(ins).fill_Grid_EWHres_MSSA_both_uptoS_RMS(end,:,:));
                        fig_name=['fill_EWHres'];
                end

                t_ss_mean(nn,:)=t_ss_mean(nn,:)+t_res;
                t_nn_mean(ss,:)=t_nn_mean(ss,:)+t_res;

                t_ss_nn(ss,nn,:)=t_res;
                t2_ss_nn_std(disp_data,ss,nn)=std(t_res);
                t2_ss_nn_rms(disp_data,ss,nn)=rms(t_res);

                s_ss_std_mean(nn,1)=s_ss_std_mean(nn,1)+sum(s_std(in2).*c11cmn12_area(in2))/sum(c11cmn12_area(in2));
                s_nn_std_mean(ss,1)=s_nn_std_mean(ss,1)+sum(s_std(in2).*c11cmn12_area(in2))/sum(c11cmn12_area(in2));
                s_ss_rms_mean(nn,1)=s_ss_rms_mean(nn,1)+sum(s_rms(in2).*c11cmn12_area(in2))/sum(c11cmn12_area(in2));
                s_nn_rms_mean(ss,1)=s_nn_rms_mean(ss,1)+sum(s_rms(in2).*c11cmn12_area(in2))/sum(c11cmn12_area(in2));

                s2_ss_nn_std(disp_data,ss,nn)=sum(s_std(in2).*c11cmn12_area(in2))/sum(c11cmn12_area(in2));
                s2_ss_nn_rms(disp_data,ss,nn)=sum(s_rms(in2).*c11cmn12_area(in2))/sum(c11cmn12_area(in2));
            end

        end

        t2_ss_mean_std(disp_data,:)=std(t_ss_mean'/numel(use_ss));
        t2_nn_mean_std(disp_data,:)=std(t_nn_mean'/numel(use_nn));
        s2_ss_mean_std(disp_data,:)=s_ss_std_mean/numel(use_ss);
        s2_nn_mean_std(disp_data,:)=s_nn_std_mean/numel(use_nn);

        t2_ss_mean_rms(disp_data,:)=rms(t_ss_mean'/numel(use_ss));
        t2_nn_mean_rms(disp_data,:)=rms(t_nn_mean'/numel(use_nn));
        s2_ss_mean_rms(disp_data,:)=s_ss_rms_mean/numel(use_ss);
        s2_nn_mean_rms(disp_data,:)=s_nn_rms_mean/numel(use_nn);
    end

    %% temporal STD (NS)

    use_ss=S_range2;
    use_nn=N_range2;

    t2_ns_std=[];t2_ns_rms=[];

    s2_ns_std=[];s2_ns_rms=[];
    for disp_data=1:5

        t_ns=zeros(numel(use_nn),numel(fill_months));

        for ii=1:numel(use_nn)

            nn=use_nn(ii);

            ss=use_ss(ii);

            %         load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 '_TP_S' num2str(ss) '_N' num2str(nn) '_M60_compare.mat']),...
            %             'A5_Centers_Compare','A5_Centers_STD_Compare');
            load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 '_TP' num2str(MNC) ...
                '_S' num2str(ss) '_N' num2str(nn) '_M' num2str(M_rec) '_sig' num2str(Noise_SigLev) ...
                '_compare.mat']),...
                'A5_Centers_Compare','A5_Centers_STD_Compare','A5_Centers_RMS_Compare');

            switch disp_data
                case 1
                    t_res=A5_Centers_Compare(ins).fill_EWH_nonoisecorrection;
                    s_std=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_uptoS_STD(end,:,:));
                    s_rms=squeeze(A5_Centers_RMS_Compare(ins).fill_Grid_EWH_MSSA_uptoS_RMS(end,:,:));
                    fig_name=['fill_EWH_nonoisecorrection'];
                case 2
                    t_res=A5_Centers_Compare(ins).fill_EWH_noisecorrection;
                    s_std=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_bothfail_uptoS_STD(end,:,:));
                    s_rms=squeeze(A5_Centers_RMS_Compare(ins).fill_Grid_EWH_MSSA_bothfail_uptoS_RMS(end,:,:));
                    fig_name=['fill_EWH_noisecorrection'];
                case 3
                    t_res=A5_Centers_Compare(ins).fill_EWH;
                    s_std=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_both_uptoS_STD(end,:,:));
                    s_rms=squeeze(A5_Centers_RMS_Compare(ins).fill_Grid_EWH_MSSA_both_uptoS_RMS(end,:,:));
                    fig_name=['fill_EWH'];
                case 4
                    t_res=A5_Centers_Compare(ins).fill_EWHsig;
                    s_std=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWHsig_MSSA_both_uptoS_STD(end,:,:));
                    s_rms=squeeze(A5_Centers_RMS_Compare(ins).fill_Grid_EWHsig_MSSA_both_uptoS_RMS(end,:,:));
                    fig_name=['fill_EWHsig'];
                case 5
                    t_res=A5_Centers_Compare(ins).fill_EWHres;
                    s_std=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWHres_MSSA_both_uptoS_STD(end,:,:));
                    s_rms=squeeze(A5_Centers_RMS_Compare(ins).fill_Grid_EWHres_MSSA_both_uptoS_RMS(end,:,:));
                    fig_name=['fill_EWHres'];
            end

            t_ns(ii,:)=t_ns(ii,:)+t_res;

            s2_ns_std(disp_data,ii)=sum(s_std(in2).*c11cmn12_area(in2))/sum(c11cmn12_area(in2));
            s2_ns_rms(disp_data,ii)=sum(s_rms(in2).*c11cmn12_area(in2))/sum(c11cmn12_area(in2));

        end

        t2_ns_std(disp_data,:)=std(t_ns');
        t2_ns_rms(disp_data,:)=rms(t_ns');

    end

    %% temporal analysis for region 1 (RMS)
    % note , the STD is very similar to RMS because the mean of GRACE is
    % subtracted already.

    %% the spatial STD of selected process
    % map_color=mymap("rainbow");
    font_Size=16;
    MarkSize=2;
    F_Position=[2100,50,1500,800];
    abc='abcdefghijklmnopqrstuvwxyz';

    use_ss=1:2:numel(S_range2);
    use_nn=1:2:numel(N_range2);

    disp_data=3;

    for ins=1

        f7=figure;

        F_gap=[.04 .03];F_marg_h=[.06 .04];F_marg_w=[.03 .08];

        ybou1=[0 12]; % SCS
        h0=tight_subplot(numel(use_ss),numel(use_nn),F_gap,F_marg_h,F_marg_w);

        % You'd better adjust this manually
        set (gcf,'Position',F_Position)

        % choose nexttile or subplot (based on your matlab version containing tiledlayout or not)
        % nexttile is recommended for testing; subplot for final publish (because of the title of colorbar maybe covered by subplots)
        %     t=tiledlayout('flow'); % nexttile method
        %     [a,b]=sizewind(N); % subplot method

        ff=0;

        for ss=use_ss

            SSF=turn_V2_four(ss,used_sta_V);
            % choose nexttile or subplot (based on your matlab version containing tiledlayout or not)
            %             nexttile % nexttile method
            %         subplot(a,b,SSF)  % subplot method

            for nn=use_nn

                load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 '_TP' num2str(MNC) ...
                    '_S' num2str(ss) '_N' num2str(nn) '_M' num2str(M_rec) '_sig' num2str(Noise_SigLev) ...
                    '_compare.mat']),'A5_Centers_STD_Compare','A5_Centers_RMS_Compare');

                ff=ff+1;

                axes(h0(ff));

                switch disp_data
                    case 1
                        r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_uptoS_STD(SSF,:,:));fig_name=[parrent_path2 '_fill_Grid_EWH_MSSA_uptoS_STD'];
                    case 2
                        r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_resid_uptoS_STD(SSF,:,:));fig_name=[parrent_path2 '_fill_Grid_EWH_MSSA_resid_uptoS_STD'];
                    case 3
                        r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_bothfail_uptoS_STD(SSF,:,:));fig_name=[parrent_path2 '_fill_Grid_EWH_MSSA_bothfail_uptoS_STD'];
                    case 4
                        r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_both_uptoS_STD(SSF,:,:));fig_name=[parrent_path2 '_fill_Grid_EWH_MSSA_both_uptoS_STD'];
                end

                %         r_res=squeeze(A5_Centers_STD(ins).fill_Grid_EWH_STD(SSF,:,:));fig_name=['fill_Grid_EWH_STD'];
                %             r_res=squeeze(A5_Centers_STD(ins).fill_Grid_EWHsig_STD(SSF,:,:));fig_name=['fill_Grid_EWHsig_STD'];
                %             r_res=squeeze(A5_Centers_STD(ins).fill_Grid_EWHres_STD(SSF,:,:));fig_name=['fill_Grid_EWHres_STD'];
                %
                %                         r_res=squeeze(A5_Centers_STD(ins).fill_Grid_EWH_uptoS_STD(SSF,:,:));fig_name=['fill_Grid_EWH_uptoS_STD'];
                %                         r_res=squeeze(A5_Centers_STD(ins).fill_Grid_EWHsig_uptoS_STD(SSF,:,:));fig_name=['fill_Grid_EWHsig_uptoS_STD'];
                %                         r_res=squeeze(A5_Centers_STD(ins).fill_Grid_EWHres_uptoS_STD(SSF,:,:));fig_name=['fill_Grid_EWHres_uptoS_STD'];
                %
                %             r_res=squeeze(A5_Centers_STD(ins).fill_Grid_EWHsig_MSSA_both_uptoS_STD(SSF,:,:));fig_name=['fill_Grid_EWHsig_MSSA_both_uptoS_STD'];
                %             r_res=squeeze(A5_Centers_STD(ins).fill_Grid_EWHres_MSSA_both_uptoS_STD(SSF,:,:));fig_name=['fill_Grid_EWHres_MSSA_both_uptoS_STD'];
                %                         r_res=squeeze(A5_Centers_STD(ins).fill_Grid_EWH_MSSA_both_uptoS_STD(SSF,:,:));fig_name=['fill_Grid_EWH_MSSA_both_uptoS_STD'];
                %                 r_res=squeeze(A5_Centers_STD(ins).fill_Grid_EWH_MSSA_bothfail_uptoS_STD(SSF,:,:));fig_name=['fill_Grid_EWH_MSSA_bothfail_uptoS_STD'];
                %             r_res=squeeze(A5_Centers_STD(ins).fill_Grid_EWH_MSSA_resid_uptoS_STD(SSF,:,:));fig_name=['fill_Grid_EWH_MSSA_resid_uptoS_STD'];

                %
                m_proj('mercator','long',[lon_c11cmn2_matl(1), lon_c11cmn2_matl(end)],'lat',[lat_c11cmn2_matl(end), lat_c11cmn2_matl(1)]);
                m_pcolor(lonlon_c11cmn2_matl,latlat_c11cmn2_matl,r_res(c11cmn_comTo2_lat,c11cmn_comTo2_lon));
                % m_coast('patch',[1 .85 .7]);
                m_coast
                %             m_grid('box','fancy','tickdir','in');

                m_line(XY2.XY_ori(:,1),XY2.XY_ori(:,2),'Linestyle','-','color','k','linewidth',1);
                hold on
                m_line(XY2.XY_buf(:,1),XY2.XY_buf(:,2),'Linestyle','--','color','r','linewidth',1);
                % xlabel(['CC(' num2str(SSF) ')'])

                %             m_text(c11cmn_combine(1)+1,c11cmn_combine(4)+5,['S=' num2str(ss) ' N=' num2str(nn) ],'color','red','fontsize',12);
                %             title(['N = ' char(use_N_name(SSF_n))],'color','red','fontsize',12)

                if nn==use_nn(1) && ss==use_ss(end)
                    m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','xtick',[90:10:130]);
                elseif nn==use_nn(1)
                    m_grid('tickdir','in','yaxislocation','left','xticklabels',[],'xtick',[90:10:130]);
                elseif ss==use_ss(end)
                    m_grid('tickdir','in','xaxislocation','bottom','yticklabels',[],'xtick',[90:10:130]);
                else
                    m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'xtick',[90:10:130])
                end

                m_text(lon_c11cmn2_matl(1)+1,lat_c11cmn2_matl(end)+4,{['S(' num2str(ss) ')=' num2str(turn_V2_four(ss,used_sta_V))], ['N(' num2str(nn) ')\approx' num2str(round(turn_MSSA2_four_ave(ss,nn)))] }, ...
                    'color','black','fontsize',font_Size-4,'verticalalignment','middle','horizontalalignment','left',...
                    'backgroundColor','white');

                %             m_line([99,99,106,106,99],[13,6,6,13,13],'color','green','linewidth',1.5);
                %             hold on

                title(['(' abc(ff) ')'], ...
                    'color','black','fontsize',font_Size,'FontName','Times New Roman')

                %             if SSF_n==1
                %                 title(use_process_name(nn),'fontsize',14)
                %             end
                %         caxis([0 10]);
                caxis(ybou1);

                set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.012,0.055]);

%                 writegeofiles(lon_combine12,lat_combine12,r_res,fullfile(output_path,[fig_name '_TP_S' num2str(ss) '_N' num2str(nn) '_SIG' num2str(Noise_SigLev) '.txt']),[],4)

            end

            %%% nexttile method
            %         t.TileSpacing='Compact';
            %         t.Padding='Compact';

            colormap(jet)
%             colormap(mymap('rainbow'))

        end

        cb = colorbar;
        cb.Position = [0.94 0.1 0.02,0.8];
        title(cb,'STD','fontsize',14);
        %%%

        whratio=abs((c11cmn(3)-c11cmn(1))/c11cmn(2)-c11cmn(4));

        tif_name7=[fig_name,'_Spatial_STD_NSchange_SIG' num2str(Noise_SigLev) '.tif'];
        print(f7,'-dtiff','-r300',['R3_' tif_name7]);

    end

    %% the spatial STD of selected process and the temporal STD for region 1 (two data combined)
    % map_color=mymap("rainbow");
    font_Size=16;
    MarkSize=2;
    F_Position=[2100,-350,1300,1100];
    abc='abcdefghijklmnopqrstuvwxyz';

    use_ss=1:2:numel(S_range2);
    use_nn=1:2:numel(N_range2);
    % use_disp=[1,3,4,5,2];
    use_disp=[1,3,2];


    disp_data=1;

    % disp_data=[1];cat_seq=1;

    ins=1;
    % for ins=1

    f7=figure;

    F_gap0=[.07 .12];F_marg_h0=[.33 .04];F_marg_w0=[.04 .36];

    ybou1=[0 25]; % SCS
    h0=tight_subplot(numel(use_ss),numel(use_nn),F_gap0,F_marg_h0,F_marg_w0);

    % You'd better adjust this manually
    set (gcf,'Position',F_Position)

    % choose nexttile or subplot (based on your matlab version containing tiledlayout or not)
    % nexttile is recommended for testing; subplot for final publish (because of the title of colorbar maybe covered by subplots)
    %     t=tiledlayout('flow'); % nexttile method
    %     [a,b]=sizewind(N); % subplot method

    ff=0;

    for ss=use_ss

        SSF=turn_V2_four(ss,used_sta_V);
        % choose nexttile or subplot (based on your matlab version containing tiledlayout or not)
        %             nexttile % nexttile method
        %         subplot(a,b,SSF)  % subplot method

        for nn=use_nn

            load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 '_TP' num2str(MNC) ...
                '_S' num2str(ss) '_N' num2str(nn) '_M' num2str(M_rec) '_sig' num2str(Noise_SigLev) ...
                '_compare.mat']),...
                'A5_Centers_STD_Compare');

            ff=ff+1;

            axes(h0(ff));

            switch disp_data(1)
                case 1
                    r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_uptoS_STD(SSF,:,:));fig_name=[parrent_path2 '_fill_Grid_EWH_MSSA_uptoS_STD'];
                case 2
                    r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_resid_uptoS_STD(SSF,:,:));fig_name=[parrent_path2 '_fill_Grid_EWH_MSSA_resid_uptoS_STD'];
                case 3
                    r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_bothfail_uptoS_STD(SSF,:,:));fig_name=[parrent_path2 '_fill_Grid_EWH_MSSA_bothfail_uptoS_STD'];
                case 4
                    r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_both_uptoS_STD(SSF,:,:));fig_name=[parrent_path2 '_fill_Grid_EWH_MSSA_both_uptoS_STD'];
            end

            %
            m_proj('mercator','long',[lon_c11cmn2_matl(1), lon_c11cmn2_matl(end)],'lat',[lat_c11cmn2_matl(end), lat_c11cmn2_matl(1)]);
            m_pcolor(lonlon_c11cmn2_matl,latlat_c11cmn2_matl,r_res(c11cmn_comTo2_lat,c11cmn_comTo2_lon));
            % m_coast('patch',[1 .85 .7]);
            m_coast
            %             m_grid('box','fancy','tickdir','in');

            m_line(XY2.XY_ori(:,1),XY2.XY_ori(:,2),'Linestyle','-','color','k','linewidth',1);
            hold on
            m_line(XY2.XY_buf(:,1),XY2.XY_buf(:,2),'Linestyle','--','color','r','linewidth',1);
            % xlabel(['CC(' num2str(SSF) ')'])

            %             m_text(c11cmn_combine(1)+1,c11cmn_combine(4)+5,['S=' num2str(ss) ' N=' num2str(nn) ],'color','red','fontsize',12);
            %             title(['N = ' char(use_N_name(SSF_n))],'color','red','fontsize',12)

            if nn==use_nn(1) && ss==use_ss(end)
                m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','xtick',[90:10:130]);
            elseif nn==use_nn(1)
                m_grid('tickdir','in','yaxislocation','left','xticklabels',[],'xtick',[90:10:130]);
            elseif ss==use_ss(end)
                m_grid('tickdir','in','xaxislocation','bottom','yticklabels',[],'xtick',[90:10:130]);
            else
                m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'xtick',[90:10:130])
            end


            %         m_text(lon_c11cmn1_matl(1)+1,lat_c11cmn1_matl(end)+5,{['\kappa_{S}:' num2str(lamda_range(ss)*100) '%'], ['\kappa_{N}:' num2str(lamda_range(nn)*100) '%'] }, ...
            %             'color','black','fontsize',font_Size-7,'verticalalignment','middle','horizontalalignment','left',...
            %             'backgroundColor','white');

            m_line([99,99,106,106,99],[13,6,6,13,13],'color',[0.99,0.99,0.99],'linewidth',1.5);
            hold on

            %         title(['(' abc(ff) ')'], ...
            %             'color','black','fontsize',font_Size,'FontName','Times New Roman')

            m_text(lon_c11cmn2_matl(end),lat_c11cmn2_matl(1)+6,['(' abc(ff) ')'],...
                'color','black','FontName','Times New Roman','fontsize',font_Size,...
                'verticalalignment','middle','horizontalalignment','center');

            %             if SSF_n==1
            %                 title(use_process_name(nn),'fontsize',14)
            %             end
            %         caxis([0 10]);
            caxis(ybou1);

            set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.012,0.055]);

            title(['EWH'], ...
                'color','black','FontName','Times New Roman','fontsize',font_Size-4)
            %             m_text(mean(lon_c11cmn2_matl),lat_c11cmn2_matl(1)-0.5,{['EWH']}, ...
            %                 'color','black','FontName','Times New Roman','fontsize',font_Size-4,'verticalalignment','middle','horizontalalignment','center',...
            %                 'backgroundColor','white');
            %         writegeofiles(lon_combine12,lat_combine12,r_res,fullfile(output_path,[fig_name '_S' num2str(lamda_range(ss)*100) '_N' num2str(lamda_range(nn)*100) '.txt']),[],4)

        end

        %         ff=ff+1;
        %%% nexttile method
        %         t.TileSpacing='Compact';
        %         t.Padding='Compact';

        colormap(jet)
%         colormap(mymap('rainbow'))

    end

    disp_data=4;

    F_gap0e=F_gap0;F_marg_h0e=F_marg_h0;F_marg_w0e=[F_marg_w0(1)+0.115 F_marg_w0(2)-0.115];

    % ybou1=[0 15]; % SCS
    h0=tight_subplot(numel(use_ss),numel(use_nn),F_gap0e,F_marg_h0e,F_marg_w0e);

    % You'd better adjust this manually
    % set (gcf,'Position',F_Position)

    % choose nexttile or subplot (based on your matlab version containing tiledlayout or not)
    % nexttile is recommended for testing; subplot for final publish (because of the title of colorbar maybe covered by subplots)
    %     t=tiledlayout('flow'); % nexttile method
    %     [a,b]=sizewind(N); % subplot method

    ff=0;

    for ss=use_ss

        SSF=turn_V2_four(ss,used_sta_V);
        % choose nexttile or subplot (based on your matlab version containing tiledlayout or not)
        %             nexttile % nexttile method
        %         subplot(a,b,SSF)  % subplot method

        for nn=use_nn

            load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 '_TP' num2str(MNC) ...
                '_S' num2str(ss) '_N' num2str(nn) '_M' num2str(M_rec) '_sig' num2str(Noise_SigLev) ...
                '_compare.mat']),...
                'A5_Centers_STD_Compare');

            ff=ff+1;

            axes(h0(ff));

            switch disp_data
                case 1
                    r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_uptoS_STD(SSF,:,:));fig_name=[parrent_path2 '_fill_Grid_EWH_MSSA_uptoS_STD'];
                case 2
                    r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_resid_uptoS_STD(SSF,:,:));fig_name=[parrent_path2 '_fill_Grid_EWH_MSSA_resid_uptoS_STD'];
                case 3
                    r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_bothfail_uptoS_STD(SSF,:,:));fig_name=[parrent_path2 '_fill_Grid_EWH_MSSA_bothfail_uptoS_STD'];
                case 4
                    r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_both_uptoS_STD(SSF,:,:));fig_name=[parrent_path2 '_fill_Grid_EWH_MSSA_both_uptoS_STD'];
            end

            %
            m_proj('mercator','long',[lon_c11cmn2_matl(1), lon_c11cmn2_matl(end)],'lat',[lat_c11cmn2_matl(end), lat_c11cmn2_matl(1)]);
            m_pcolor(lonlon_c11cmn2_matl,latlat_c11cmn2_matl,r_res(c11cmn_comTo2_lat,c11cmn_comTo2_lon));
            % m_coast('patch',[1 .85 .7]);
            m_coast
            %             m_grid('box','fancy','tickdir','in');

            m_line(XY2.XY_ori(:,1),XY2.XY_ori(:,2),'Linestyle','-','color','k','linewidth',1);
            hold on
            m_line(XY2.XY_buf(:,1),XY2.XY_buf(:,2),'Linestyle','--','color','r','linewidth',1);
            % xlabel(['CC(' num2str(SSF) ')'])

            %             m_text(c11cmn_combine(1)+1,c11cmn_combine(4)+5,['S=' num2str(ss) ' N=' num2str(nn) ],'color','red','fontsize',12);
            %             title(['N = ' char(use_N_name(SSF_n))],'color','red','fontsize',12)

            if ss==use_ss(end)
                m_grid('tickdir','in','xaxislocation','bottom','yticklabels',[],'xtick',[90:10:130]);
            else
                m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'xtick',[90:10:130])
            end

            %         m_text(lon_c11cmn1_matl(1)+1,lat_c11cmn1_matl(end)+5,{['S' num2str(ss)], ['N' num2str(nn)] }, ...
            %             'color','black','fontsize',font_Size-7,'verticalalignment','middle','horizontalalignment','left',...
            %             'backgroundColor','white');

            m_line([99,99,106,106,99],[13,6,6,13,13],'color',[0.99,0.99,0.99],'linewidth',1.5);
            hold on

            %         title(['(' abc(ff) ')'], ...
            %             'color','black','fontsize',font_Size,'FontName','Times New Roman')

            %             if SSF_n==1
            %                 title(use_process_name(nn),'fontsize',14)
            %             end
            %         caxis([0 10]);
            caxis(ybou1);

            set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.012,0.055]);

            title(['EWH (noise removal)'], ...
                'color','black','FontName','Times New Roman','fontsize',font_Size-4);
            %             m_text(mean(lon_c11cmn2_matl),lat_c11cmn2_matl(1)-0.5,{['EWH (noise removal)']}, ...
            %                 'color','black','FontName','Times New Roman','fontsize',font_Size-4,'verticalalignment','middle','horizontalalignment','center',...
            %                 'backgroundColor','white');
            m_text(lon_c11cmn2_matl(1)-8,lat_c11cmn2_matl(end)+4.5,{['S(' num2str(ss) ')=' num2str(turn_V2_four(ss,used_sta_V))], ['N(' num2str(nn) ')\approx' num2str(round(turn_MSSA2_four_ave(ss,nn)))] }, ...
                'color','black','fontsize',font_Size-5,'verticalalignment','middle','horizontalalignment','left',...
                'backgroundColor','white','FontName','Times New Roman');


            %         writegeofiles(lon_combine12,lat_combine12,r_res,fullfile(output_path,[fig_name '_S' num2str(lamda_range(ss)*100) '_N' num2str(lamda_range(nn)*100) '.txt']),[],4)

        end

        %         ff=ff+1;
        %%% nexttile method
        %         t.TileSpacing='Compact';
        %         t.Padding='Compact';

                    colormap(jet)
%         colormap(mymap('rainbow'))

    end

    cb = colorbar;
    cb.Location = 'South';
    cb.Position = [0.10 0.275 0.6,0.015];
    title(cb,'RMS','fontsize',14);
    %%%

    whratio=abs((c11cmn(3)-c11cmn(1))/c11cmn(2)-c11cmn(4));

    colorr_blind=[0,114,178;204,121,167;0,158,115;230,159,0;213,94,0;240,228,66;86,180,233]/255; % Blind-friendly colors

    ybou1=[0 10]; % Yangtze
    ybou2=[0 4.5]; % Yangtze
    xbou=[1-0.5 S_range2(end)+0.5];

    % use_disp=[1,3,4,5,2];
    use_disp=[1,3,2];

    % notice that here length(use_ss) should be equal to length(use_nn)
    xticklabelS={};xticklabelN={};xticklabelSN={};
    xticklabelS_gap={};xticklabelN_gap={};
    for i=1:length(use_ss)
        xticklabelS{i}=['S(' num2str(use_ss(i)) ')'];
        xticklabelS_gap{use_ss(i)}=['S(' num2str(use_ss(i)) ')'];
    end
    for i=1:length(use_nn)
        xticklabelN{i}=['N(' num2str(use_nn(i)) ')'];
        xticklabelN_gap{use_nn(i)}=['N(' num2str(use_nn(i)) ')'];
    end
    for i=1:length(use_nn)
        xticklabelSN{i}=['S(' num2str(use_ss(i)) '), N(' num2str(use_nn(i)) ')'];
    end

    % F1-F3
    Max_ratio_nn=zeros(length(use_nn),2);
    % use_ss_nn_rms=t2_ss_nn_rms;use_sn_rms=t2_ns_rms;ybou1=[0 5];RMS_T='GRMS';
    use_ss_nn_rms=s2_ss_nn_rms;use_sn_rms=s2_ns_rms;ybou1=[0 20];RMS_T='WRMS';

    disp_sn=3; %signal rms
    use_ss_nn_square=squeeze(use_ss_nn_rms(disp_sn,:,:));


    F_gap1=[F_gap0(1) 0.05];F_marg_h1=[.06 .77];F_marg_w1=[F_marg_w0(1)+0.02 F_marg_w0e(2)+0.02];

    h1=tight_subplot(1,length(use_nn),F_gap1,F_marg_h1,F_marg_w1);

    xbou1=[0.5,MNC+0.5];

    for i=1:length(use_nn)

        ff=ff+1;

        num_pic=i;

        axes(h1(num_pic));


        for disp=use_disp
            p1(disp)=plot(S_range2,use_ss_nn_rms(disp,:,use_nn(i)),'Marker','none','LineStyle','-',...
                'LineWidth',2,'color',colorr_blind(disp,:));
            hold on
        end

        % xlabel('\lambda (%)')
        if i==1
            %             ylabel('$\overline{\mathrm{RMS}}$','Interpreter','latex')
            ylabel(RMS_T)
        else
            %         yticklabels([]);
        end

        t=text(S_range2(1),ybou1(2)*4/5,['N=N(' num2str(use_nn(i)) ')'],'FontName','Times New Roman','FontSize',font_Size-3);
        t.EdgeColor='k';
        t.LineStyle='-';
        t.LineWidth=1;

        ylim(ybou1);
        xlim(xbou1);
        xticks(use_ss)
        xticklabels(xticklabelS)

        title(['(' abc(ff) ')'], ...
            'color','black','fontsize',font_Size,'FontName','Times New Roman')

        set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.012,0.055]);

    end

    F_gap2=[.06 .03];F_marg_h2=[.32 .03];F_marg_w2=[.81 .02];

    h2=tight_subplot(length(use_ss),1,F_gap2,F_marg_h2,F_marg_w2);

    Max_ratio_ss=zeros(length(use_ss),2);
    for i=1:length(use_ss)

        ff=ff+1;
        num_pic=i;

        axes(h2(num_pic));

        pn=[];
        for disp=use_disp
            p1(disp)=plot(N_range2,reshape(use_ss_nn_rms(disp,use_ss(i),:),1,[]),'Marker','none','LineStyle','-',...
                'LineWidth',2,'color',colorr_blind(disp,:));
            hold on
            pn=[pn p1(disp)];
        end

        ylabel(RMS_T)

        t=text(N_range2(1),ybou1(2)*4/5,['S=S(' num2str(use_ss(i)) ')'],'FontName','Times New Roman','FontSize',font_Size-3);
        t.EdgeColor='k';
        t.LineStyle='-';
        t.LineWidth=1;

        ylim(ybou1);
        xlim(xbou1);
        xticks(use_nn)
        xticklabels(xticklabelN)

        title(['(' abc(ff) ')'], ...
            'color','black','fontsize',font_Size,'FontName','Times New Roman')

        set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.012,0.055]);
    end

    F_gap3=F_gap0;F_marg_h3=F_marg_h1;F_marg_w3=[F_marg_w2(1)-0.01 F_marg_w2(2)];

    h2=tight_subplot(1,1,F_gap3,F_marg_h3,F_marg_w3);

    axes(h2(1));

    ff=ff+1;

    % F3 (new)
    SHM=SHeatmap(use_ss_nn_square,'Format','sq');
    SHM=SHM.draw();
    SHM.setText('FontSize',font_Size+2);

    for i=1:size(use_ss_nn_square,1)
        for j=1:size(use_ss_nn_square,2)
            if use_ss_nn_square(i,j)==max(use_ss_nn_square,[],'all')
                SHM.setTextMN(i,j,'String','x')
            else
                SHM.setTextMN(i,j,'String','')
            end
            %             if use_ss_nn_square(i,j)==max(diag(use_ss_nn_square),[],'all')
            %                 SHM.setTextMN(i,j,'String','o')
            %             end
        end
    end

    title(['(' abc(ff) ')'], ...
        'color','black','fontsize',font_Size,'FontName','Times New Roman')
    % SHM.setText('FontSize',font_Size+2);
%     set(gca,'colormap',mymap('coolwarm'));
    set(gca,'colormap','jet');

    ax=gca;
    ax.XTickLabel=xticklabelN_gap;
    ax.YTickLabel=xticklabelS_gap;
    ax.FontSize=14;

    cb=colorbar;
    cb.Location='eastoutside';
    cb.Position = [0.965 0.06 0.012 0.17];
    %     title(cb,{'$\overline{\mathrm{RMS}}$'},'fontsize',14,'Interpreter','latex');
    title(cb,{RMS_T},'fontsize',14,'Interpreter','latex');

    set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.012,0.055]);

    legend_text={'EWH','EWH (noise)','EWH (noise removal)','Fitted signal of EWH (noise removal)','Unfitted signal of EWH (noise removal)'};
    l1=legend(pn,legend_text(use_disp),...
        'box','off','Location',[0.15,0.013,0.7,0.02],'FontSize',font_Size,'NumColumns',3);

    tif_name7=[fig_name,'_Spatial_Temporal_twoSTD_NSchange_SIG' num2str(Noise_SigLev) '.tif'];
    print(f7,'-dtiff','-r500',['R3_' tif_name7]);

    Square_Seq{SIG}=use_ss_nn_square;
end

%%
Square_Seq2=Square_Seq;

save([Code_Version '_' suffix_area2 '_' suffix_XY_use2  '_S' num2str(Max_S) suffix_bou  '_Data.mat'],'Square_Seq2',...
    'SIG_use','use_ss','use_nn','MNC','M_rec','RMS_T','c11cmn','c11cmn_combine')
%%
load([Code_Version '_' suffix_area2 '_' suffix_XY_use2  '_S' num2str(Max_S) suffix_bou  '_Data.mat'])

%%
% notice that here length(use_ss) should be equal to length(use_nn)
xticklabelS={};xticklabelN={};xticklabelSN={};
xticklabelS_gap={};xticklabelN_gap={};
for i=1:length(use_ss)
    xticklabelS{i}=['S(' num2str(use_ss(i)) ')'];
    xticklabelS_gap{use_ss(i)}=['S(' num2str(use_ss(i)) ')'];
end
for i=1:length(use_nn)
    xticklabelN{i}=['N(' num2str(use_nn(i)) ')'];
    xticklabelN_gap{use_nn(i)}=['N(' num2str(use_nn(i)) ')'];
end
for i=1:length(use_nn)
    xticklabelSN{i}=['S(' num2str(use_ss(i)) '), N(' num2str(use_nn(i)) ')'];
end

f8=figure;

font_Size=16;
MarkSize=2;
abc='abcdefghijklmnopqrstuvwxyz';

F_Position=[2100,50,1100,300];
F_gap0=[.04 .03];F_marg_h0=[.1 .1];F_marg_w0=[.05 .1];

ybou1=[0 12]; % SCS
h0=tight_subplot(1,length(SIG_use),F_gap0,F_marg_h0,F_marg_w0);

% You'd better adjust this manually
set (gcf,'Position',F_Position)

ff=0;
for jj=1:length(SIG_use)
    axes(h0(jj));

    ff=ff+1;

    use_ss_nn_square_SIG=Square_Seq2{jj};
    % F3 (new)
    SHM=SHeatmap(use_ss_nn_square_SIG,'Format','sq');
    SHM=SHM.draw();
    SHM.setText('FontSize',font_Size+2);

    for i=1:size(use_ss_nn_square_SIG,1)
        for j=1:size(use_ss_nn_square_SIG,2)
            if use_ss_nn_square_SIG(i,j)==max(use_ss_nn_square_SIG,[],'all')
                SHM.setTextMN(i,j,'String','x')
                SIG_ss(jj)=i;
                SIG_nn(jj)=j;
            else
                SHM.setTextMN(i,j,'String','')
            end
            if use_ss_nn_square_SIG(i,j)==max(diag(use_ss_nn_square_SIG),[],'all')
                SHM.setTextMN(i,j,'String','o')
            end
        end
    end

    title(['(' abc(ff) ') ' num2str(100-SIG_use(jj)*100) '%'], ...
        'color','black','fontsize',font_Size+4,'FontName','Times New Roman')
    % SHM.setText('FontSize',font_Size+2);
%     set(gca,'colormap',mymap('coolwarm'));
    set(gca,'colormap','jet');

    ax=gca;
    ax.XTickLabel=xticklabelN_gap;
    if jj==1
        ax.YTickLabel=xticklabelS_gap;
    else
        ax.YTickLabel={};
    end
    ax.FontSize=14;

    clim(ybou1)

    if jj==length(SIG_use)
        cb=colorbar;
        cb.Location='eastoutside';
        cb.Position = [0.94 0.15 0.015 0.7];
        %         title(cb,{'$\overline{\mathrm{RMS}}$'},'fontsize',14,'Interpreter','latex');
        title(cb,RMS_T,'fontsize',14);
    else
        colorbar off
    end

    set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.012,0.055]);

end

% tif_name8=[fig_name,'_Heatmap_4SIG.tif'];
% print(f8,'-dtiff','-r300',['R3_' tif_name8]);

%
fprintf('Then we determine the optimal N and S.');
% SIG_ss(1:length(SIG_use))=5;
% SIG_nn(1:length(SIG_use))=5;

save([Code_Version '_' suffix_area2 '_' suffix_XY_use2  '_S' num2str(Max_S) suffix_bou '_Data'],'SIG_ss','SIG_nn','-append');

%%
f9=figure;

disp_SIG=[1,2,3];

F_Position=[2100,50,1100,500];

% You'd better adjust this manually
set (gcf,'Position',F_Position)

% tif_name8=[fig_name,'_Heatmap_4SIG.tif'];
% print(f8,'-dtiff','-r300',['R3_' tif_name8]);

%
col_gap=[0 205 205]/255;col_est=[153 204 255]/255;col_out=[128 0 128]/255;
col_art=[238 173 14]/255;
% map_color=mymap("rainbow");
font_Size=16;
MarkSize=2;
% F_Position=[2100,-250,1200,500];
F_gap0=[.04 .03];F_marg_h0=[.2 .05];F_marg_w0=[.07 .04];
abc='abcdefghijklmnopqrstuvwxyz';

% box plot for the energy distribution
colorr_blind=[213,94,0;204,121,167;0,158,115;230,159,0;213,94,0;240,228,66;0,114,178;86,180,233]/255; % Blind-friendly colors

% f9=figure;

left_color=[0 0 0];
right_color=[255 0 0]/255;
set(f9,'defaultAxesColorOrder',[left_color;right_color])

ybou1=[0 1.5];
h0=tight_subplot(1,length(disp_SIG),F_gap0,F_marg_h0,F_marg_w0);


ff=0;
for jj=disp_SIG

    load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 '_TP' num2str(MNC) ...
        '_S' num2str(SIG_ss(jj)) '_N' num2str(SIG_nn(jj)) '_M' num2str(M_rec) '_sig' num2str(SIG_use(jj)) ...
        '_compare.mat']),...
        'A5_Centers_Compare','S');

    ff=ff+1;

    axes(h0(ff));

    used_eva=A5_Centers_Compare(1).fill_MSSAcoffs_eva;

    %         yyaxis left
    x=[];
    % notice that we do not want IB to show here.
    %         if S_ol>S
    %             for i=1:S
    %                 x{i}=num2str(i);
    %             end
    %             x{S+1}='IB';
    %             X = categorical(x);
    %             X = reordercats(X,x);
    %         else
    %             X=1:S_ol;
    %
    %         end

    X=1:S;
    Y=[used_eva(1:S,[3,4,6]), used_eva(1:S,[8,9,11]) , ones(size(used_eva(1:S,1))) - used_eva(1:S,1)];

    bb=bar(X,Y,1,'stacked');
    set(bb(1),'FaceColor',colorr_blind(1,:));
    set(bb(2),'FaceColor',colorr_blind(3,:));
    % set(bb(3),'FaceColor',colorr_blind(2,:));
    set(bb(3),'FaceColor',colorr_blind(4,:));

    set(bb(4),'FaceColor',colorr_blind(1,:),'FaceAlpha',0.3);
    set(bb(5),'FaceColor',colorr_blind(3,:),'FaceAlpha',0.3);
    % set(bb(3),'FaceColor',colorr_blind(2,:));
    set(bb(6),'FaceColor',colorr_blind(4,:),'FaceAlpha',0.3);

    %         hatchfill2(bb(4),'single','HatchAngle',60,'hatchcolor',colorr_blind(1,:),'HatchLineWidth',1.5);
    %         hatchfill2(bb(5),'single','HatchAngle',60,'hatchcolor',colorr_blind(3,:),'HatchLineWidth',1.5);
    %         hatchfill2(bb(6),'single','HatchAngle',60,'hatchcolor',colorr_blind(4,:),'HatchLineWidth',1.5);
    %         for b = 4:6
    %             %     bb(b).FaceColor = [0.8,0.8,0.8];
    %             bb(b).FaceColor = 'none';
    %         end


    % set(bb(4),'FaceColor',colorr_blind(1,:));
    % set(bb(5),'FaceColor',colorr_blind(3,:));
    % % set(bb(3),'FaceColor',colorr_blind(2,:));
    % set(bb(6),'FaceColor',colorr_blind(4,:));

    set(bb(7),'FaceColor',[169,169,169]/255);

    ylim(ybou1);
    xlim([0,S+1])

    grid on

    if ff==1
        ylabel({['sum of normalized eigenvalues'], ['in M-SSA ($\sum\widetilde{\lambda}$)']},'Interpreter','latex')
        %         legend([bb(1),bb(2),bb(3),bb(7)],{'long-term','annual','short-term','residual'}, ...
        %             'NumColumns',3,'Location','northeast','box','off')
        legend([bb(1),bb(2),bb(3),bb(4),bb(5),bb(6),bb(7)],{'long-term','annual',...
            'short-term','long-term (noise)','annual (noise)',...
            'short-term (noise)','residual'}, ...
            'NumColumns',7,'Position',[0.2,0.02,0.6,0.05],'box','off')

        %             legend([bb(1),bb(2),bb(3)],{'long-term','annual',...
        %                 'short-term','long-term (noise)','annual (noise)',...
        %                 'short-term (noise)','residual'}, ...
        %                 'NumColumns',2,'Location','north','box','off')
    elseif ff==2
        %             legend([bb(4),bb(5),bb(6)],{'long-term (noise)','annual (noise)',...
        %                 'short-term (noise)','residual'}, ...
        %                 'NumColumns',2,'Location','north','box','off')
    else
        %             legend([bb(7)],{'residual'}, ...
        %                 'NumColumns',2,'Location','north','box','off')
    end
    xlabel('Spherical slepian coefficients')

    if ff>1
        yticklabels([]);
    end

    title({['\gamma = ' num2str((1-SIG_use(jj))*100) '%']},...
        'FontSize',font_Size-2,'color','r')

    set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.004,0.035]);

    xticks(get(gca, 'xtick'))
    set(gca, 'xticklabel', get(gca, 'xtick'))


end

F_gap0=[.04 .2];F_marg_h0=[.71 .08];F_marg_w0=[.22 .066];

ybou1=[0 15]; % SCS
h0=tight_subplot(1,length(disp_SIG),F_gap0,F_marg_h0,F_marg_w0);
ff=0;
for jj=disp_SIG

    ff=ff+1;

    axes(h0(ff));

    use_ss_nn_square_SIG=Square_Seq2{jj};
    % F3 (new)
    SHM=SHeatmap(use_ss_nn_square_SIG,'Format','sq');
    SHM=SHM.draw();
    SHM.setText('FontSize',font_Size-2);

    for i=1:size(use_ss_nn_square_SIG,1)
        for j=1:size(use_ss_nn_square_SIG,2)
            if use_ss_nn_square_SIG(i,j)==max(use_ss_nn_square_SIG,[],'all')
                SHM.setTextMN(i,j,'String','x')
            else
                SHM.setTextMN(i,j,'String','')
            end
            %             if use_ss_nn_square_SIG(i,j)==max(diag(use_ss_nn_square_SIG),[],'all')
            %                 SHM.setTextMN(i,j,'String','o')
            %             end
        end
    end

    %     title(['(' abc(ff) ') ' num2str(100-SIG_use(jj)*100) '%'], ...
    %         'color','black','fontsize',font_Size+4,'FontName','Times New Roman')
    % SHM.setText('FontSize',font_Size+2);
    set(gca,'colormap',mymap('coolwarm'));
    set(gca,'colormap','jet');

    ax=gca;
    %     ax.XTickLabel=xticklabelN_gap;
    %     ax.YTickLabel=xticklabelS_gap;
    %     if jj==1
    %         ax.YTickLabel=xticklabelS_gap;
    %     else
    %         ax.YTickLabel={};
    %     end
    ax.FontSize=font_Size-6;

    clim(ybou1)

    %     if ff==length(disp_SIG)
    cb=colorbar;
    cb.Location='eastoutside';
    %         cb.Position = [0.92 0.73 0.015 0.17];
    %         title(cb,{'$\overline{\mathrm{RMS}}$'},'fontsize',14,'Interpreter','latex');
    title(cb,{'WRMS'},'fontsize',font_Size-6,'Interpreter','latex');
    %     else
    %         colorbar off
    %     end

    %     set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.012,0.055]);

end

%%
% % here is a bug fix
% load(fullfile(getenv('IFILES'),'Bug1_usedate.mat'),'common_use_months')
% use_member=ismember(A5_Centers_Compare(1).fill_months,common_use_months);
% missing_months=A5_Centers_Compare(1).fill_months(~use_member);

Noise_SigLev=0.05;
R95=load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 '_TP' num2str(MNC) ...
    '_S' num2str(SIG_ss(1)) '_N' num2str(SIG_nn(1)) '_M' num2str(M_rec) '_sig' num2str(Noise_SigLev) ...
    '_compare.mat']),'A5_Centers_Compare');
R95.A5_Centers_Compare.name='CJ0323_R95';

Noise_SigLev=0.1;
% load(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' Attach_ALL '_S' num2str(lamda_range(ss)*100) '_N' num2str(lamda_range(nn)*100) '_compare.mat']))
R90=load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 '_TP' num2str(MNC) ...
    '_S' num2str(SIG_ss(2)) '_N' num2str(SIG_nn(2)) '_M' num2str(M_rec) '_sig' num2str(Noise_SigLev) ...
    '_compare.mat']),'A5_Centers_Compare');
R90.A5_Centers_Compare.name='CJ0323_R90';

Noise_SigLev=0.3;
% load(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' Attach_ALL '_S' num2str(lamda_range(ss)*100) '_N' num2str(lamda_range(nn)*100) '_compare.mat']))
R70=load(fullfile(getenv('IFILES'),['Results_' parrent_path2],['Main_' Attach_ALL2 '_TP' num2str(MNC) ...
    '_S' num2str(SIG_ss(3)) '_N' num2str(SIG_nn(3)) '_M' num2str(M_rec) '_sig' num2str(Noise_SigLev) ...
    '_compare.mat']),'A5_Centers_Compare');
R70.A5_Centers_Compare.name='CJ0323_R70';

A5_Centers_Compare=R95.A5_Centers_Compare;
A5_Centers_Compare(2)=R90.A5_Centers_Compare;
A5_Centers_Compare(3)=R70.A5_Centers_Compare;
% A5_Centers_Compare(4)=R99.A5_Centers_Compare;

%% Description of variables in A5_Centers_Compare
%name   |   The name of the product
%c11cmn |   The [lon,lat] boundary for presenting the results of study of interest and surrounding regions (should include the entire study region of interest)
%lon    |   The longitude
%lat    |   The latitude
%data_year_beg  |   The begining year of the data period
%process        |   The post-process for GRACE data (e.g., 'RL06' means Rlease 06 of GRACE; '60' means the maxiumum degree of spherical harmonics; 'land' indicates the study region is land; 'Pelt17' means the GIA models)
%data_date      |   The date for each data. Because we use multiple data centers, we do not provide a consistent data date.
%fill_date      |   The date with gap-filling
%fill_month     |   The month with gap-filling
%use_date       |   The date used for interpolate the gap (i.e., the commonly used data across different GRACE/GRACE-FO data centers).
%use_month      |   The month used for interpolate the gap 
%missing_date   |   The date missing
%missing_month  |   The month missing
%fill_MSSAcoffs_eva         |   Sum of Eigenvalues of MSSA for different cases (reference to 'Lite_MSSA_V1.m')
%fill_MSSAcoffs_RCTest      |   Information for noise detection for each spherical Slepian coefficietns (SSC) during the MSSA.
%fill_EWH_Slepian           |   The area-weighted time series of equivalent water height (EWH) (unit: cm) only after spherical Slepian method
%fill_EWH_nonoisecorrection |   The area-weighted time series of EWH only after the cutoff process of MSSA (without noise process)
%fill_EWH_noisecorrection   |   The area-weighted time series of EWH recognized as noise during the noise detection process of MSSA
%fill_EWH                   |   The area-weighted time series of EWH after the cutoff and noise correction of MSSA
%fill_EWHsig                |   The area-weighted signals of time series 'fill_EWH' extracted by eight-coefficient least-squares
%fill_EWHres                |   The area-weighted residuals of time series 'fill_EWH' extracted by eight-coefficient least-squares
%fill_uptoS_EWH_nonoisecorrection   |   The area-weighted time series of 'fill_EWH_nonoisecorrection' but accumulating by the increase of truncation number S
%fill_uptoS_EWH_noisecorrection     |   The area-weighted time series of 'fill_EWH_noisecorrection' but accumulating by the increase of truncation number S
%fill_uptoS_EWH                     |   The area-weighted time series of 'fill_EWH' but accumulating by the increase of truncation number S
%fill_uptoS_EWHsig                  |   The area-weighted time series of 'fill_EWHsig' but accumulating by the increase of truncation number S
%fill_uptoS_EWHres                  |   The area-weighted time series of 'fill_EWHres' but accumulating by the increase of truncation number S
%fill_Grid_EWH_Slepian              |   The spatial maps of EWH only after spherical Slepian method
%fill_Grid_EWH_nonoisecorrection    |   The spatial maps of EWH only after the cutoff process of MSSA (without noise process)
%fill_Grid_EWH                      |   The spatial maps of EWH after the cutoff and noise correction of MSSA
%fill_Grid_EWHsig                   |   The spatial maps of signals after the cutoff and noise correction of MSSA extracted by eight-coefficient least-squares
%fill_Grid_EWHres                   |   The spatial maps of residuals after the cutoff and noise correction of MSSA extracted by eight-coefficient least-squares
%fill_MASS                          |   The area-weighted time series of mass (unit: Gt) after the cutoff and noise correction of MSSA
%fill_MASSsig                       |   The area-weighted signals of time series 'fill_MASS' extracted by eight-coefficient least-squares
%fill_MASSres                       |   The area-weighted residuals of time series 'fill_MASS' extracted by eight-coefficient least-squares
%fill_Grid_MASS                     |   The spatial maps of mass (unit: Gt) after the cutoff and noise correction of MSSA
%fill_Grid_MASSsig                  |   The spatial maps of signals of 'fill_Grid_MASS' extracted by eight-coefficient least-squares
%fill_Grid_MASSres                  |   The spatial maps of residuals of 'fill_Grid_MASS' extracted by eight-coefficient least-squares
%fill_MSL_Slepian           |   The area-weighted time series of non-steric sea level anomalies (NSLA) (unit: cm) only after spherical Slepian method
%fill_MSL                   |   The area-weighted time series of NSLA after the cutoff and noise correction of MSSA
%fill_MSLsig                |   The area-weighted signals of time series 'fill_NSLA' extracted by eight-coefficient least-squares
%fill_MSLres                |   The area-weighted residuals of time series 'fill_NSLA' extracted by eight-coefficient least-squares
%fill_Grid_MSL_Slepian              |   The spatial maps of MSL only after spherical Slepian method
%fill_Grid_MSL                      |   The spatial maps of MSL after the cutoff and noise correction of MSSA
%fill_Grid_MSLsig                   |   The spatial maps of signals after the cutoff and noise correction of MSSA extracted by eight-coefficient least-squares
%fill_Grid_MSLres                   |   The spatial maps of residuals after the cutoff and noise correction of MSSA extracted by eight-coefficient least-squares
%fill_IBdeltacoeffs                 |   The inverted barometer correction (reference to 'Lite_MSSA_V1.m')
%fill_GADfilldelta                  |   The spatial maps of GAD (reference to 'Lite_MSSA_V1.m')


%%

%%
fitwhat=[3 365.25 365.25/2]; %linear trend and an annual and a semiannual harmonic term

% for ins=[1:11,13]
%     A5_Centers_Compare(ins).tt_fit_beg=1;A5_Centers_Compare(ins).tt_fit_end=252; % we want to delete one year in the begining and end
% end
% A5_Centers_Compare(12).tt_fit_beg=1;A5_Centers_Compare(12).tt_fit_end=237; % we want to delete one year in the begining and end

for ins=[1:3]
    A5_Centers_Compare(ins).tt_fit_beg=13;A5_Centers_Compare(ins).tt_fit_end=240; % we want to delete one year in the begining and end
end

fitwhat=[3 365.25 365.25/2]; %linear trend and an annual and a semiannual harmonic term

%temporal
for ins=1:3

    %     if ins==3
    %         A5_Centers_Compare(ins).tt_fit_beg=13;A5_Centers_Compare(ins).tt_fit_end=237; % we want to delete one year in the begining and end (for Gauer)
    %     else
    %         A5_Centers_Compare(ins).tt_fit_beg=13;A5_Centers_Compare(ins).tt_fit_end=240; % we want to delete one year in the begining and end
    %     end

    fit_data_length=A5_Centers_Compare(ins).tt_fit_beg:A5_Centers_Compare(ins).tt_fit_end;

    fill_dates=A5_Centers_Compare(ins).fill_dates(fit_data_length);
    fill_MSL=A5_Centers_Compare(ins).fill_MSL(fit_data_length);
    fill_MSL_resid=A5_Centers_Compare(ins).fill_MSLres(fit_data_length);

    A5_Centers_Compare(ins).fill_MSL_Std=std(fill_MSL);

    % kk=kk-mean(kk);

    [signal_kk,resid_kk,~,~,...
        Amplitude,Phase]=periodic_analysis_m(fill_MSL',...
        fill_dates,fitwhat); % sin(wt-fi)
    Phase=pi/2-Phase; % be consistent with previous studies for cos(wt+fi)
    Phase(Phase<0)=Phase(Phase<0)+2*pi; % -180~-180 to 0~360
    A5_Centers_Compare(ins).fill_MSL_Amp=Amplitude;
    A5_Centers_Compare(ins).fill_MSL_Pha=Phase/pi*180;

    A5_Centers_Compare(ins).fill_MSL_signal_Std=std(signal_kk);
    A5_Centers_Compare(ins).fill_MSL_resid_Std=std(resid_kk);

    [fit_MSL,delta_MSL,totalparams_MSL,paramerrors_MSL] = timeseriesfit([fill_dates' fill_MSL'], ...
        var(resid_kk),1,1);
    A5_Centers_Compare(ins).fill_MSL_slope=totalparams_MSL(2)*365.25;
    A5_Centers_Compare(ins).fill_MSL_fit = [fill_dates' fit_MSL delta_MSL];
    A5_Centers_Compare(ins).fill_MSL_paramerrors = paramerrors_MSL(2)*365.25; % this (365.25) is special for 'timeseriesfit' results.

    % seasonal and deseasonal (trend is removed)
    fill_MSL_detrend=fill_MSL-fit_MSL';
    fill_MSL_deseason=fill_MSL_detrend;fill_MSL_season=zeros(1,12);
    for j=1:12
        fill_MSL_season(j)=mean(fill_MSL_detrend(j:12:end));

        fill_MSL_deseason(j:12:end)=fill_MSL_deseason(j:12:end)-fill_MSL_season(j);
    end
    A5_Centers_Compare(ins).fill_MSL_season=fill_MSL_season;
    A5_Centers_Compare(ins).fill_MSL_deseason=fill_MSL_deseason;
end

%% spatial

for ins=1:3

    fit_data_length=A5_Centers_Compare(ins).tt_fit_beg:A5_Centers_Compare(ins).tt_fit_end;

    fill_dates=A5_Centers_Compare(ins).fill_dates(fit_data_length);
    fill_MSL=A5_Centers_Compare(ins).fill_Grid_MSL(fit_data_length,:,:);
    %     fill_MSL_resid=A5_Centers_Compare(ins).fill_Grid_MSLres;

    [t,m,n]=size(fill_MSL);

    A5_Centers_Compare(ins).fill_Grid_NaN=nan(m,n);
    A5_Centers_Compare(ins).fill_MSL_Std=nan(m,n);
    A5_Centers_Compare(ins).fill_Grid_MSL_Amp=nan(numel(fitwhat)-1,m,n);
    A5_Centers_Compare(ins).fill_Grid_MSL_Pha=nan(numel(fitwhat)-1,m,n);
    A5_Centers_Compare(ins).fill_Grid_MSLsig=nan(t,m,n);
    A5_Centers_Compare(ins).fill_Grid_MSLres=nan(t,m,n);
    A5_Centers_Compare(ins).fill_Grid_MSL_signal_Std=nan(m,n);
    A5_Centers_Compare(ins).fill_Grid_MSL_resid_Std=nan(m,n);
    A5_Centers_Compare(ins).fill_Grid_MSL_slope=nan(m,n);
    A5_Centers_Compare(ins).fill_Grid_MSL_fit ={};
    A5_Centers_Compare(ins).fill_Grid_MSL_paramerrors = nan(m,n) ; % this (365.25) is special for 'timeseriesfit' results.

    A5_Centers_Compare(ins).fill_Grid_MSL_season=nan(12,m,n);
    A5_Centers_Compare(ins).fill_Grid_MSL_deseason=nan(t,m,n);
    A5_Centers_Compare(ins).fill_Grid_MSL_4seasons=nan(4,m,n);

    A5_Centers_Compare(ins).fill_Grid_MSLsig_season=nan(12,m,n);
    A5_Centers_Compare(ins).fill_Grid_MSLsig_deseason=nan(t,m,n);
    A5_Centers_Compare(ins).fill_Grid_MSLsig_4seasons=nan(4,m,n);

    A5_Centers_Compare(ins).fill_Grid_MSLres_season=nan(12,m,n);
    A5_Centers_Compare(ins).fill_Grid_MSLres_deseason=nan(t,m,n);
    A5_Centers_Compare(ins).fill_Grid_MSLres_4seasons=nan(4,m,n);

    for mm=1:m

        mm

        for nn=1:n

            fill_one_MSL=fill_MSL(:,mm,nn)';

            if isnan(fill_one_MSL)
                A5_Centers_Compare(ins).fill_Grid_NaN(mm,nn)=sum(isnan(fill_one_MSL));
                continue
            end

            A5_Centers_Compare(ins).fill_MSL_Std(mm,nn)=std(fill_one_MSL);

            % kk=kk-mean(kk);

            [signal_kk,resid_kk,~,~,...
                Amplitude,Phase]=periodic_analysis_m(fill_one_MSL',...
                fill_dates,fitwhat); % sin(wt-fi)
            Phase=pi/2-Phase; % be consistent with previous studies for cos(wt+fi)
            Phase(Phase<0)=Phase(Phase<0)+2*pi; % -180~-180 to 0~360
            A5_Centers_Compare(ins).fill_Grid_MSL_Amp(:,mm,nn)=Amplitude;
            A5_Centers_Compare(ins).fill_Grid_MSL_Pha(:,mm,nn)=Phase;

            A5_Centers_Compare(ins).fill_Grid_MSLsig(:,mm,nn)=signal_kk;
            A5_Centers_Compare(ins).fill_Grid_MSLres(:,mm,nn)=resid_kk;

            A5_Centers_Compare(ins).fill_Grid_MSL_signal_Std(mm,nn)=std(signal_kk);
            A5_Centers_Compare(ins).fill_Grid_MSL_resid_Std(mm,nn)=std(resid_kk);

            [fit_MSL,delta_MSL,totalparams_MSL,paramerrors_MSL] = timeseriesfit([fill_dates' fill_one_MSL'], ...
                var(resid_kk),1,1);
            A5_Centers_Compare(ins).fill_Grid_MSL_slope(mm,nn)=totalparams_MSL(2)*365.25;
            A5_Centers_Compare(ins).fill_Grid_MSL_fit{mm,nn} = [fill_dates' fit_MSL delta_MSL];
            A5_Centers_Compare(ins).fill_Grid_MSL_paramerrors(mm,nn) = paramerrors_MSL(2)*365.25; % this (365.25) is special for 'timeseriesfit' results.

            A5_Centers_Compare(ins).fill_Grid_MSL_trend(:,mm,nn)=fit_MSL;

            % seasonal
            fill_MSL_detrend=fill_one_MSL-fit_MSL';
            A5_Centers_Compare(ins).fill_Grid_MSL_detrend(:,mm,nn)=fill_MSL_detrend;

            fill_MSL_deseason=fill_MSL_detrend;fill_MSL_season=zeros(1,12);
            for j=1:12
                fill_MSL_season(j)=mean(fill_MSL_detrend(j:12:end));

                fill_MSL_deseason(j:12:end)=fill_MSL_deseason(j:12:end)-fill_MSL_season(j);
            end
            A5_Centers_Compare(ins).fill_Grid_MSL_season(:,mm,nn)=fill_MSL_season;
            A5_Centers_Compare(ins).fill_Grid_MSL_deseason(:,mm,nn)=fill_MSL_deseason;

            A5_Centers_Compare(ins).fill_Grid_MSL_4seasons(1,mm,nn)=mean(fill_MSL_season(3:5));
            A5_Centers_Compare(ins).fill_Grid_MSL_4seasons(2,mm,nn)=mean(fill_MSL_season(6:8));
            A5_Centers_Compare(ins).fill_Grid_MSL_4seasons(3,mm,nn)=mean(fill_MSL_season(9:11));
            A5_Centers_Compare(ins).fill_Grid_MSL_4seasons(4,mm,nn)=mean(fill_MSL_season([12,1:2]));

            % seasonal signal
            fill_MSLsig_deseason=signal_kk;fill_MSLsig_season=zeros(1,12);
            fill_MSLres_deseason=resid_kk;fill_MSLres_season=zeros(1,12);
            for j=1:12
                fill_MSLsig_season(j)=mean(signal_kk(j:12:end));
                fill_MSLres_season(j)=mean(resid_kk(j:12:end));

                fill_MSLsig_deseason(j:12:end)=fill_MSLsig_deseason(j:12:end)-fill_MSLsig_season(j);
                fill_MSLres_deseason(j:12:end)=fill_MSLres_deseason(j:12:end)-fill_MSLres_season(j);
            end
            A5_Centers_Compare(ins).fill_Grid_MSLsig_season(:,mm,nn)=fill_MSLsig_season;
            A5_Centers_Compare(ins).fill_Grid_MSLsig_deseason(:,mm,nn)=fill_MSLsig_deseason;
            A5_Centers_Compare(ins).fill_Grid_MSLres_season(:,mm,nn)=fill_MSLres_season;
            A5_Centers_Compare(ins).fill_Grid_MSLres_deseason(:,mm,nn)=fill_MSLres_deseason;

            A5_Centers_Compare(ins).fill_Grid_MSLsig_4seasons(1,mm,nn)=mean(fill_MSLsig_season(3:5));
            A5_Centers_Compare(ins).fill_Grid_MSLsig_4seasons(2,mm,nn)=mean(fill_MSLsig_season(6:8));
            A5_Centers_Compare(ins).fill_Grid_MSLsig_4seasons(3,mm,nn)=mean(fill_MSLsig_season(9:11));
            A5_Centers_Compare(ins).fill_Grid_MSLsig_4seasons(4,mm,nn)=mean(fill_MSLsig_season([12,1:2]));

            A5_Centers_Compare(ins).fill_Grid_MSLres_4seasons(1,mm,nn)=mean(fill_MSLres_season(3:5));
            A5_Centers_Compare(ins).fill_Grid_MSLres_4seasons(2,mm,nn)=mean(fill_MSLres_season(6:8));
            A5_Centers_Compare(ins).fill_Grid_MSLres_4seasons(3,mm,nn)=mean(fill_MSLres_season(9:11));
            A5_Centers_Compare(ins).fill_Grid_MSLres_4seasons(4,mm,nn)=mean(fill_MSLres_season([12,1:2]));

        end
    end

end

%%
% f10=figure;
% Temporal (three figures)
ax_name=["STPC(95%)","STPC(90%)","STPC(70%)"];

for ins=[1:3]
    A5_Centers_Compare(ins).tt_fit_beg=13;A5_Centers_Compare(ins).tt_fit_end=240; % we want to delete one year in the begining and end
end

comp_data=[1:3];
% comp_data=[6,8,9,10,5];
% comp_data=[1:13];
tt_fit_year=2004:2022;

fitwhat=[3 365.25 365.25/2]; %linear trend and an annual and a semiannual harmonic term

tt_use=floor(A5_Centers_Compare(1).use_months/12)+A5_Centers_Compare(1).data_year_beg+mod(A5_Centers_Compare(1).use_months,12)/12-1/24;
% tt_fil=floor(fill_months/12)+data_year_beg+mod(fill_months,12)/12-1/24;
tt_lea=floor(A5_Centers_Compare(1).missing_months/12)+A5_Centers_Compare(1).data_year_beg+mod(A5_Centers_Compare(1).missing_months,12)/12-1/24;

col_gap=[0 205 205]/255;col_est=[153 204 255]/255;col_out=[128 0 128]/255;
col_art=[238 173 14]/255;
colorr_blind=[0,114,178;204,121,167;0,158,115;230,159,0;213,94,0;240,228,66;86,180,233]/255; % Blind-friendly colors

colorr_use=colorr_blind;

f1=figure
% plot the basis
F_gap=[.04 .03];F_marg_h=[.06 .04];F_marg_w=[.05 .03];

font_Size=16;
MarkSize=2;
F_Position=[100,50,1500,600];
abc='abcdefghijklmnopqrstuvwxyz';

ybou1=[-15,15];
h0=tight_subplot(1,1,F_gap,F_marg_h,F_marg_w);
axes(h0(1));

for j=1:numel(tt_lea)
    leak_tt_neibour=[tt_lea(j)-1/24,tt_lea(j),tt_lea(j)+1/24];
    %         leak_MSL_neibour=[Total_MSL_EST_reconst(missing_months(j)-1),...
    %             Total_MSL_EST_reconst(missing_months(j)), ...
    %             Total_MSL_EST_reconst(missing_months(j)+1)];
    %         p2=plot(leak_tt_neibour,leak_MSL_neibour,'r','linewidth',2);

    he_gap=area(leak_tt_neibour([1,3]),[ybou1(2);ybou1(2)],ybou1(1));
    he_gap.FaceColor = col_gap;
    he_gap.EdgeColor = col_gap;
    he_gap.FaceAlpha = 0.4;
    he_gap.EdgeAlpha = 0.7;

    hold on
end

pn=[];lg={};
ff=1;
for ins=comp_data

    fit_data_length=A5_Centers_Compare(ins).tt_fit_beg:A5_Centers_Compare(ins).tt_fit_end;
    fill_dates=A5_Centers_Compare(ins).fill_dates(fit_data_length);
    fill_MSL=A5_Centers_Compare(ins).fill_MSL(fit_data_length);

    tt_fil=floor(A5_Centers_Compare(ins).fill_months/12)+A5_Centers_Compare(ins).data_year_beg+mod(A5_Centers_Compare(ins).fill_months,12)/12-1/24;
    tt_use=floor(A5_Centers_Compare(ins).fill_months(fit_data_length)/12)+A5_Centers_Compare(ins).data_year_beg+...
        mod(A5_Centers_Compare(ins).fill_months(fit_data_length),12)/12-1/24;

    p1(ff)=plot(tt_use,A5_Centers_Compare(ins).fill_MSL(fit_data_length),'linewidth',2,'color',colorr_use(ins,:));
    hold on
    p2(ff)=plot(tt_use,A5_Centers_Compare(ins).fill_MSL_fit(:,2),'color',colorr_use(ins,:),'linewidth',1);
    hold on
    %     p3(ff)=plot(tt_fil,(A5_Centers_Compare(ins).fill_MSL_fit(:,2)+A5_Centers_Compare(ins).fill_MSL_fit(:,3)),'m','linestyle','--','linewidth',0.7);
    %     hold on
    %     p4(ff)=plot(tt_fil,(A5_Centers_Compare(ins).fill_MSL_fit(:,2)-A5_Centers_Compare(ins).fill_MSL_fit(:,3)),'m','linestyle','--','linewidth',0.7);
    %     hold on
    %
    pn=[pn p1(ff)];
    lg{ff}=[char(ax_name(ins)) ': ' num2str(A5_Centers_Compare(ins).fill_MSL_slope,'%.3f') ' \pm ' ...
        num2str(A5_Centers_Compare(ins).fill_MSL_paramerrors,'%.3f') ' cm/yr'] ;

    ff=ff+1;
end

xticks(tt_fit_year)
xlim([tt_fit_year(1)-0.5 tt_fit_year(end)+1.5])
ylim(ybou1)

ylabel('Area-weighted MSL (cm)')
xlabel('Time (year)')

ld=legend(pn,lg,'box','off','Location','Northwest','FontSize',font_Size-1,'NumColumns',3);

set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.004,0.035]);

set(gcf,'Position',F_Position)

tif_name1=['Temporal_MSL_compare_other_SIG.tif'];
print(f1,'-dtiff','-r300',fullfile(fig_path_ALL2,tif_name1));
