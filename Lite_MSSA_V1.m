% This code is to reconstruct the mass and EWH (MSL) time series based on
% Multiple Singular Spectrum Analysis (M-SSA) for spherical Slepian
% Functions, suitable for advance users because there are some paramters
% for adjustment if necessary.
% Code version: 1 (2nd round RSE revision).
% Code Note:
%   Main function:
%       (1) Fit Gaps through the M-SSA;
%       (2) Reconstruct the Spherical slepian coefficients through M-SSA based
%           both the K-S test and Lilliefors test.
%       (3) Provide results based on different p values
%   Main features:
%       (1) do the trend fit based on reconstructed slepian coefficients
%       (2*) adaptive for all regions
%
% Tested on 9.13.0.2049777 (R2022b)
% Last modified by zhongtian.ma@connect.polyu.hk, 12/05/2025

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

% boundary restriction of S or N
threshold_S=[S_bou, 1-S_bou]; % restriction on extremely high or low lambda
threshold_N=[N_bou, 1-N_bou];

% redo some initial data? (default 0)
redo=0;

% for different regions, the buffer zone you should decided reference to simulation

%%%% (case 1) information

parrent_path=Code_Version;
c11cmn_combine=c11cmn; 

vv=char(Dataproduct(2));
suffix_area=[vv(3:4) '_' Area];

suffix_XY_use=[num2str(Lwindow) '_' buffer_str '_' num2str(Radius)];
%%%%

%  (notice that this is fixed to be SSF=50)
suffix_str=[suffix_area '_' suffix_XY_use '_S' num2str(Max_S)];

if S_bou==0 && N_bou==0
    suffix_bou=[];
elseif S_bou==0
    suffix_bou=['_Hn' num2str(N_bou)];
elseif N_bou==0
    suffix_bou=['_Hs' num2str(S_bou)];
else
    suffix_bou=['_Hs' num2str(S_bou) 'n' num2str(N_bou) ];
end

% please rewrite the name of this path accordingly to the output of 'V5r1.m'
Attach_ALL=[char(use_institu(intit_num+1)) '_' suffix_str   suffix_bou];
txt_path_ALL=fullfile(getenv('IFILES'),['Text_' parrent_path],Attach_ALL);
fig_path_ALL=fullfile(getenv('IFILES'),['Figure_' parrent_path],Attach_ALL);
for ins=1:intit_num
    Attach_each(ins)=string([char(use_institu(ins)) '_' suffix_str]);
    txt_path_each(ins)=string(fullfile(getenv('IFILES'),['Text_' parrent_path],[char(Attach_each(ins)) suffix_bou]));
    fig_path_each(ins)=string(fullfile(getenv('IFILES'),['Figure_' parrent_path],[char(Attach_each(ins))  suffix_bou]));

    if ~exist(char(txt_path_each(ins)),'dir')
        mkdir(char(txt_path_each(ins)))
    end

    if ~exist(char(fig_path_each(ins)),'dir')
        mkdir(char(fig_path_each(ins)))
    end

end
Attach_each(intit_num+1)=string(Attach_ALL);
txt_path_each(intit_num+1)=string(txt_path_ALL);
fig_path_each(intit_num+1)=string(fig_path_ALL);
%%%%

XY_path=fullfile(getenv('IFILES'),['Results_' parrent_path],['XY_' ...
    char(use_institu(1)) '_' suffix_area '_' suffix_XY_use]);
load(XY_path)

if ~exist(txt_path_ALL,'dir')
    mkdir(txt_path_ALL)
end

if ~exist(fig_path_ALL,'dir')
    mkdir(fig_path_ALL)
end


%% Load the corresponding data for each institution

if redo || ~exist(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' Attach_ALL 'c11.mat']),'file')

    for ins=1:intit_num

        dataname=fullfile(getenv('IFILES'),['Results_' parrent_path],['MainData_' char(Attach_each(ins)) '.mat'])

        load(dataname)

        for i=1:length(Slept_dates)
            monthstart = datenum([data_year_beg...
                i 1]); % assume begin from January
            monthend = datenum([data_year_beg...
                i+1 1]);
            model_dates(i)=(monthstart + monthend)/2;
        end
        %% Boundary output
        fid=fopen(fullfile(char(txt_path_each(ins)),'boundary_XY.txt'),'wt');
        fprintf(fid,'>\n');
        for i=1:size(XY,1)
            fprintf(fid,'%f\t%f\n',XY(i,1),XY(i,2));
        end
        fclose(fid);
        % fprintf(fid,'%f\t%f\n',XY(1,1),XY(1,2));

        figure
        plot(XY(:,1),XY(:,2))
        hold on
        plot(XY_ori(:,1),XY_ori(:,2))

        close all

        S_shannon=round((Lwindow+1)^2*spharea(XY_buf));
        %% Calcualte GRACE

        %     save(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' FIG_Attach '.mat']),'Table1','Table1_un')

        Dataprocess=Dataproduct;
        Dataprocess(1)=[];

        A5_Centers_Sort(ins).name=char(use_institu(ins));
        A5_Centers_Sort(ins).data_year_beg=data_year_beg;A5_Centers_Sort(ins).process=Dataprocess;
        A5_Centers_Sort(ins).data_dates=data_dates;A5_Centers_Sort(ins).fill_dates=fill_dates;
        A5_Centers_Sort(ins).use_dates=use_dates;A5_Centers_Sort(ins).use_months=use_months';
        A5_Centers_Sort(ins).missing_dates=missing_dates;A5_Centers_Sort(ins).missing_months=missing_months;

        A5_Centers_Sort(ins).V=V;

        % basin-all
        A5_Centers_Sort(ins).fill_deltacoffs=fill_deltacoffs;
        A5_Centers_Sort(ins).fill_residdelta=fill_residdelta;A5_Centers_Sort(ins).fill_signaldelta=fill_signaldelta;
        A5_Centers_Sort(ins).fill_MASS=fill_MASS; % Gt
        A5_Centers_Sort(ins).fill_MASSsig=fill_MASSsig; % Gt
        A5_Centers_Sort(ins).fill_MASSres=fill_MASSres; % Gt

        % SSF coefficients
        A5_Centers_Sort(ins).use_residcoffs=use_residcoffs;
        A5_Centers_Sort(ins).use_sleptdelta=use_sleptdelta;

        A5_Centers_Fit(ins).use_MASSfit=use_MASSfit; % Gt
        A5_Centers_Fit(ins).data_MASSslope=data_MASSparams(2)*365.25; % Gt
        A5_Centers_Fit(ins).data_MASSparamerrors = data_MASSparamerrors(2); % Gt

        % this is the correct version of old A5 one.
%         S=N;S_sig=N_sig; % we want to use S to replace the old 'N' to denote the number of SSF
        if strcmp(Dataproduct{4},'ocean')

            S_ol=S+1;

            A5_Centers_Sort(ins).fill_IBdeltacoffs=fill_IBresiddelta+fill_IBsignaldelta; % mmH2O (~kg/m2)
            A5_Centers_Sort(ins).fill_IBresiddelta=fill_IBresiddelta; % mmH2O
            A5_Centers_Sort(ins).fill_IBsignaldelta=fill_IBsignaldelta; % mmH2O

        else

            S_ol=S;

        end

        A5_Centers_Sort(ins).fill_Grid_EWH=fill_Grid_EWH; % cm
        A5_Centers_Sort(ins).fill_Grid_EWHsig=fill_Grid_EWHsig; % cm
        A5_Centers_Sort(ins).fill_Grid_EWHres=fill_Grid_EWHres; % cm

        %     for SHC (actually it is not the main purpose here)
        disp('Calculate the SHC, it may take a while. please wait.')
        [potcoffs,~,~]=grace2plmt_m(Dataproduct{1},Dataproduct{2},Dataproduct{3},'SD', ...
            0,Dataproduct{4},Dataproduct{5});

        SHCpotcoffsdelta=potcoffs(1:size(potcoffs,1),:,:) - [repmat(mean(potcoffs(1:size(potcoffs,1),:,:),1),size(potcoffs,1),1)];

        SHCpotcoffsdelta(:,:,1:2)=potcoffs(:,:,1:2);

        [SHC_time,SHC_coefficient_SHC_values]=size(SHCpotcoffsdelta);
        sample_SHC=squeeze(SHCpotcoffsdelta(1,:,:));

        % reoder for first fit
        SHC_reoderfit=[];
        for i=1:SHC_time
            SHC_reoderfit(i,:)=reshape(squeeze(SHCpotcoffsdelta(i,:,3:4)),[],1);
        end

        fitwhat=[3 365.25 365.25/2]; %linear trend and an annual and a semiannual harmonic term

        % only use the fit function in this code to get the residual of IB
        [data_SHCsignalcoffs,data_SHCresidcoffs,data_SHCftests,data_SHCextravalues]=slept2resid_m(SHC_reoderfit,...
            Slept_dates,fitwhat); % mm (values equivalent to kg/m^2)
        use_SHCsignalcoffs=data_SHCsignalcoffs(data2use_monIndex,:);
        use_SHCresidcoffs=data_SHCresidcoffs(data2use_monIndex,:);

        fill_nmonths=length(fill_dates);use_nmonths = length(use_months);
        if any(missing_months)
            % complete estiamted signal with filled leakage values
            fill_SHCsignalcoffs_com=zeros(length(fill_dates),size(use_SHCsignalcoffs,2));
            fill_SHCresidcoffs_com=zeros(length(fill_dates),size(use_SHCresidcoffs,2));
            for i=1:size(use_SHCsignalcoffs,1)
                fill_SHCsignalcoffs_com(use_months(i),:)=use_SHCsignalcoffs(i,:);
                fill_SHCresidcoffs_com(use_months(i),:)=use_SHCresidcoffs(i,:);
            end

            for i=1:size(missing_months,2)
                fill_SHCsignalcoffs_com(missing_months(i),:)=data_SHCextravalues(i,:);
                %             SHCresid_com(leakage(i),:)=extravalues_resid(i,:);
                fill_SHCresidcoffs_com(missing_months(i),:)=0;
            end
            fill_SHCsignaldelta=fill_SHCsignalcoffs_com(1:fill_nmonths,:) - repmat(mean(use_SHCsignalcoffs(1:use_nmonths,:),1),fill_nmonths,1);
            fill_SHCresiddelta=fill_SHCresidcoffs_com(1:fill_nmonths,:) - repmat(mean(use_SHCresidcoffs(1:use_nmonths,:),1),fill_nmonths,1);

        else
            fill_SHCsignaldelta=use_SHCsignalcoffs(1:use_nmonths,:) - repmat(mean(use_SHCsignalcoffs(1:use_nmonths,:),1),use_nmonths,1);
            fill_SHCresiddelta=use_SHCresidcoffs(1:use_nmonths,:) - repmat(mean(use_SHCresidcoffs(1:use_nmonths,:),1),use_nmonths,1);
        end
        fill_SHCfilldelta=fill_SHCsignaldelta+fill_SHCresiddelta; % unit kg/m^2

        A5_Centers_Sort(ins).fill_SHCfilldelta=fill_SHCfilldelta;
        A5_Centers_Sort(ins).fill_SHCresiddelta=fill_SHCresiddelta;
        A5_Centers_Sort(ins).fill_SHCsignaldelta=fill_SHCsignaldelta;

        if strcmp(Dataproduct{4},'ocean')
            %     for GAD SHC (actually it is not the main purpose here)
            % this calculation is for subtracting from Gauer et al.,
            disp('Calculate the SHC for GAD, it may take a while. please wait.')
            [potcoffs_gsm,~,~]=grace2plmt_m(Dataproduct{1},Dataproduct{2},Dataproduct{3},'SD', ...
                0,'land',Dataproduct{5});

            % recover the GAD by removing the so-called EWH(GSM) from OBP
            % (GSM+GAD)
            SHCpotcoffsdelta_gad=potcoffs(1:size(potcoffs,1),:,:) - potcoffs_gsm(1:size(potcoffs_gsm,1),:,:) - ...
                [repmat(mean(potcoffs(1:size(potcoffs,1),:,:) - potcoffs_gsm(1:size(potcoffs_gsm,1),:,:),1),size(potcoffs_gsm,1),1)];

            SHCpotcoffsdelta_gad(:,:,1:2)=potcoffs_gsm(:,:,1:2);

            [SHC_time,SHC_coefficient_SHC_values]=size(SHCpotcoffsdelta_gad);
            sample_SHC=squeeze(SHCpotcoffsdelta_gad(1,:,:));

            [r_GAD_Sample,lon,lat,Plm,degres]=plm2xyz(sample_SHC,1,c11cmn,60);

            %             figure
            %             imagesc(r_GAD_Sample)

            % reoder for first fit
            SHC_reoderfit=[];
            for i=1:SHC_time
                SHC_reoderfit(i,:)=reshape(squeeze(SHCpotcoffsdelta_gad(i,:,3:4)),[],1);
            end

            fitwhat=[3 365.25 365.25/2]; %linear trend and an annual and a semiannual harmonic term

            % only use the fit function in this code to get the residual of IB
            [data_SHCsignalcoffs,data_SHCresidcoffs,data_SHCftests,data_SHCextravalues]=slept2resid_m(SHC_reoderfit,...
                Slept_dates,fitwhat); % mm (values equivalent to kg/m^2)
            use_SHCsignalcoffs=data_SHCsignalcoffs(data2use_monIndex,:);
            use_SHCresidcoffs=data_SHCresidcoffs(data2use_monIndex,:);

            fill_nmonths=length(fill_dates);use_nmonths = length(use_months);
            if any(missing_months)
                % complete estiamted signal with filled leakage values
                fill_SHCsignalcoffs_com=zeros(length(fill_dates),size(use_SHCsignalcoffs,2));
                fill_SHCresidcoffs_com=zeros(length(fill_dates),size(use_SHCresidcoffs,2));
                for i=1:size(use_SHCsignalcoffs,1)
                    fill_SHCsignalcoffs_com(use_months(i),:)=use_SHCsignalcoffs(i,:);
                    fill_SHCresidcoffs_com(use_months(i),:)=use_SHCresidcoffs(i,:);
                end

                for i=1:size(missing_months,2)
                    fill_SHCsignalcoffs_com(missing_months(i),:)=data_SHCextravalues(i,:);
                    %             SHCresid_com(leakage(i),:)=extravalues_resid(i,:);
                    fill_SHCresidcoffs_com(missing_months(i),:)=0;
                end
                fill_SHCsignaldelta=fill_SHCsignalcoffs_com(1:fill_nmonths,:) - repmat(mean(use_SHCsignalcoffs(1:use_nmonths,:),1),fill_nmonths,1);
                fill_SHCresiddelta=fill_SHCresidcoffs_com(1:fill_nmonths,:) - repmat(mean(use_SHCresidcoffs(1:use_nmonths,:),1),fill_nmonths,1);

            else
                fill_SHCsignaldelta=use_SHCsignalcoffs(1:use_nmonths,:) - repmat(mean(use_SHCsignalcoffs(1:use_nmonths,:),1),use_nmonths,1);
                fill_SHCresiddelta=use_SHCresidcoffs(1:use_nmonths,:) - repmat(mean(use_SHCresidcoffs(1:use_nmonths,:),1),use_nmonths,1);
            end
            fill_SHCfilldelta=fill_SHCsignaldelta+fill_SHCresiddelta; % unit kg/m^2

            A5_Centers_Sort(ins).fill_GADfilldelta=fill_SHCfilldelta;
            A5_Centers_Sort(ins).fill_GADresiddelta=fill_SHCresiddelta;
            A5_Centers_Sort(ins).fill_GADsignaldelta=fill_SHCsignaldelta;
        end

    end

    %
    % search for the common months with real GRACE measurements among products
    common_use_months=1:numel(fill_dates);
    for ins=1:intit_num
        common_use_months_index=ismember(A5_Centers_Sort(ins).use_months,common_use_months);
        common_use_months=A5_Centers_Sort(ins).use_months(common_use_months_index);
    end

    %
    A5_Centers_Sort(intit_num+1).name=char(use_institu(intit_num+1));
    A5_Centers_Sort(intit_num+1).data_year_beg=data_year_beg;A5_Centers_Sort(intit_num+1).process=Dataprocess;
    A5_Centers_Sort(intit_num+1).data_dates=[]; % default using CSR (actually this is wrong)
    A5_Centers_Sort(intit_num+1).fill_dates=fill_dates;
    A5_Centers_Sort(intit_num+1).use_dates=fill_dates(common_use_months); % default using common months % bug 1 fix
    A5_Centers_Sort(intit_num+1).use_months=common_use_months;  % default using common months % bug 1 fix
    A5_Centers_Sort(intit_num+1).missing_dates=fill_dates(setdiff(fill_months,fill_months(common_use_months))); % default using common months
    A5_Centers_Sort(intit_num+1).missing_months=fill_months(setdiff(fill_months,fill_months(common_use_months)))'; % default using common months


    fill_deltacoffs=0;fill_residdelta=0;fill_signaldelta=0;
    fill_MASS=0;fill_MASSsig=0;fill_MASSres=0;
    use_sleptdelta=0;use_residcoffs=0;

    fill_SHCfilldelta=0;fill_SHCresiddelta=0;fill_SHCsignaldelta=0;

    fill_Grid_EWH=0;fill_Grid_EWHres=0;fill_Grid_EWHsig=0;
    for ins=1:intit_num
        fill_deltacoffs=fill_deltacoffs+A5_Centers_Sort(ins).fill_deltacoffs;
        fill_residdelta=fill_residdelta+A5_Centers_Sort(ins).fill_residdelta;
        fill_signaldelta=fill_signaldelta+A5_Centers_Sort(ins).fill_signaldelta;

        fill_MASS=fill_MASS+A5_Centers_Sort(ins).fill_MASS;
        fill_MASSsig=fill_MASSsig+A5_Centers_Sort(ins).fill_MASSsig;
        fill_MASSres=fill_MASSres+A5_Centers_Sort(ins).fill_MASSres;

        fill_Grid_EWH=fill_Grid_EWH+A5_Centers_Sort(ins).fill_Grid_EWH;
        fill_Grid_EWHres=fill_Grid_EWHres+A5_Centers_Sort(ins).fill_Grid_EWHres;
        fill_Grid_EWHsig=fill_Grid_EWHsig+A5_Centers_Sort(ins).fill_Grid_EWHsig;

        fill_SHCfilldelta=fill_SHCfilldelta+A5_Centers_Sort(ins).fill_SHCfilldelta;
        fill_SHCresiddelta=fill_SHCresiddelta+A5_Centers_Sort(ins).fill_SHCresiddelta;
        fill_SHCsignaldelta=fill_SHCsignaldelta+A5_Centers_Sort(ins).fill_SHCsignaldelta;

        % transform the use month to be consistent with the first instution
        % (common months here)
        common_index=ismember(A5_Centers_Sort(ins).use_months, common_use_months);

        use_sleptdelta=use_sleptdelta+A5_Centers_Sort(ins).use_sleptdelta(common_index,:);
        use_residcoffs=use_residcoffs+A5_Centers_Sort(ins).use_residcoffs(common_index,:);
    end
    A5_Centers_Sort(intit_num+1).fill_deltacoffs=fill_deltacoffs/intit_num;
    A5_Centers_Sort(intit_num+1).fill_residdelta=fill_residdelta/intit_num;
    A5_Centers_Sort(intit_num+1).fill_signaldelta=fill_signaldelta/intit_num;

    A5_Centers_Sort(intit_num+1).fill_MASS=fill_MASS/intit_num; % Gt
    A5_Centers_Sort(intit_num+1).fill_MASSsig=fill_MASSsig/intit_num; % Gt
    A5_Centers_Sort(intit_num+1).fill_MASSres=fill_MASSres/intit_num; % Gt
    A5_Centers_Sort(intit_num+1).V=V;

    A5_Centers_Sort(intit_num+1).fill_Grid_EWH=fill_Grid_EWH/intit_num; % cm
    A5_Centers_Sort(intit_num+1).fill_Grid_EWHres=fill_Grid_EWHres/intit_num; % cm
    A5_Centers_Sort(intit_num+1).fill_Grid_EWHsig=fill_Grid_EWHsig/intit_num; % cm

    A5_Centers_Sort(intit_num+1).fill_SHCfilldelta=fill_SHCfilldelta/intit_num;
    A5_Centers_Sort(intit_num+1).fill_SHCresiddelta=fill_SHCresiddelta/intit_num;
    A5_Centers_Sort(intit_num+1).fill_SHCsignaldelta=fill_SHCsignaldelta/intit_num;

    A5_Centers_Sort(intit_num+1).use_sleptdelta=use_sleptdelta/intit_num;
    A5_Centers_Sort(intit_num+1).use_residcoffs=use_residcoffs/intit_num;

    if strcmp(Dataproduct{4},'ocean')
        fill_IBdeltacoffs=0;fill_IBresiddelta=0;fill_IBsignaldelta=0;
        for ins=1:intit_num
            fill_IBdeltacoffs=fill_IBdeltacoffs+A5_Centers_Sort(ins).fill_IBdeltacoffs;
            fill_IBresiddelta=fill_IBresiddelta+A5_Centers_Sort(ins).fill_IBresiddelta;
            fill_IBsignaldelta=fill_IBsignaldelta+A5_Centers_Sort(ins).fill_IBsignaldelta;
        end
        A5_Centers_Sort(intit_num+1).fill_IBdeltacoffs=fill_IBdeltacoffs/intit_num;
        A5_Centers_Sort(intit_num+1).fill_IBresiddelta=fill_IBresiddelta/intit_num;
        A5_Centers_Sort(intit_num+1).fill_IBsignaldelta=fill_IBsignaldelta/intit_num;

        fill_GADfilldelta=0;fill_GADresiddelta=0;fill_GADsignaldelta=0;
        for ins=1:intit_num
            fill_GADfilldelta=fill_GADfilldelta+A5_Centers_Sort(ins).fill_GADfilldelta;
            fill_GADresiddelta=fill_GADresiddelta+A5_Centers_Sort(ins).fill_GADresiddelta;
            fill_GADsignaldelta=fill_GADsignaldelta+A5_Centers_Sort(ins).fill_GADsignaldelta;
        end
        A5_Centers_Sort(intit_num+1).fill_GADfilldelta=fill_GADfilldelta/intit_num;
        A5_Centers_Sort(intit_num+1).fill_GADresiddelta=fill_GADresiddelta/intit_num;
        A5_Centers_Sort(intit_num+1).fill_GADsignaldelta=fill_GADsignaldelta/intit_num;
    end

    % calculate the mean trend
    [Cab_ave] = slepresid2cov(A5_Centers_Sort(intit_num+1).use_residcoffs);
    alphavarall_ave = functionintegrals*Cab_ave(1:S,1:S)*functionintegrals';

    [use_MASSave_CC,~,~,~,~,~]=integral_fit(S,A5_Centers_Sort(intit_num+1).use_dates,functionintegrals, ...
        A5_Centers_Sort(intit_num+1).use_sleptdelta,Cab_ave);

    use_MASSave=sum(use_MASSave_CC);

    [use_ave_fit,use_ave_delta,use_ave_params,use_ave_paramerrors] = timeseriesfit([A5_Centers_Sort(intit_num+1).use_dates' use_MASSave'],...
        alphavarall_ave,1,1);
    % Make a matrix for the line, and 95% confidence in the fit
    use_MASSavefit = [A5_Centers_Sort(intit_num+1).use_dates' use_ave_fit use_ave_delta];
    data_MASSaveparamerrors = use_ave_paramerrors*365.25;

    A5_Centers_Fit(intit_num+1).use_MASSfit=use_MASSavefit; % Gt
    A5_Centers_Fit(intit_num+1).data_MASSslope=use_ave_params(2)*365.25;
    A5_Centers_Fit(intit_num+1).data_MASSparamerrors = data_MASSaveparamerrors(2);

    % degree and order of SHC
    SHC_degord=SHCpotcoffsdelta(1,:,1:2);

    save('Bug1_usedate','common_use_months');

    save(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' Attach_ALL 'c11.mat']), ...
        'A5_Centers_Sort', 'A5_Centers_Fit','c11cmn','Lwindow','S','S_sig','S_ol','S_shannon','CC','SHC_degord','BasinArea')

else

    % directly load this variables

    load(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' Attach_ALL 'c11.mat']));

end

%%
data_year_beg=A5_Centers_Sort(1).data_year_beg;
fill_dates=A5_Centers_Sort(1).fill_dates;fill_months=1:numel(fill_dates);

V=A5_Centers_Sort(1).V;

% S_range=sum(V>(1-lamda_range));

% time wasting
[r_example,lon,lat,Plm,degres]=plm2xyz(CC{1},1,c11cmn_combine,60);
%     c11cmn2=[0.5 89.5 359.5 -89.5];
%     [r_example2,lon,lat,Plm,degres]=plm2xyz(CC{1},1,c11cmn2,60);
%         c11cmn3=[0.5 89.5 60.5 -49.5];
%     [r_example3,lon,lat,Plm,degres]=plm2xyz(CC{1},1,c11cmn3,60);
[lonlon,latlat]=meshgrid(lon,lat);

if length(S_sig)>1
    error('We do not want to use Gaussian Smooth here.')
    S1=S_sig(1);S2=S_sig(2);
    for j=1:S1
        [r_record(j,:,:),lon,lat,Plm,degres]=plm2xyz(CC{j},1,c11cmn_combine,60);
    end
    for j=S1+1:S2
        CC_smo=plm_Gausssmooth(CC{j},Radius);
        [r_record(j,:,:),lon,lat,Plm,degres]=plm2xyz(CC_smo,1,c11cmn_combine,60);
    end

else

    for j=1:S
        [r_record(j,:,:),lon,lat,Plm,degres]=plm2xyz(CC{j},1,c11cmn_combine,60);
    end
end

% CC power
CC_power=zeros(Lwindow,Max_S);
for nn=1:Max_S
    CC_N=CC{nn};

    CC_beg=1;
    for CCl=1:Lwindow
        CC_beg=CC_beg+CCl;
        for CCm=0:CCl
            CC_power(CCl,nn)=CC_power(CCl,nn)+CC_N(CC_beg+CCm,3)^2+CC_N(CC_beg+CCm,4)^2;
        end
    end
end

%%
% F_gap=[.03 .03];F_marg_h=[.02 .02];F_marg_w=[.03 .02];
% F_Position=[2100,-250,1300,1200];
%
% figure
% ha=tight_subplot(5,10,F_gap,F_marg_h,F_marg_w);
%
% % You'd better adjust this manually
% set (gcf,'Position',F_Position)
%
% for i=1:Max_S
%     axes(ha(i));
%
%     plot(1:60,CC_power(:,i),'LineWidth',1.5)
%     title(['SSF = ' num2str(i)],'FontSize',14,'FontName','Times New Roman')
% end

%% turning point of V (Slepian part)

disp('Now we determine different S by searching for breaking points of the distribution of lambda.')

% so we choose 'mean' as the chosen one
option_sta=["mean","rms","std","linear"];

% use_S_TP=find(V>0.1 & V<0.9);
use_S_TP=find(V>threshold_S(1) & V<threshold_S(2));

used_sta_V=4; % what kind of turning point do you want? (default is RSS)
MaxNumChanges_V=Turning_number; % how many turning point do you want?

turn_V_four=nan(MaxNumChanges_V,4);
for i=1:4
    used_sta=char(option_sta(i));
    % whole SCS
    [turn_V] = findchangepts(V(use_S_TP), 'Statistic', used_sta, 'MaxNumChanges', MaxNumChanges_V);

    leng_V=length(turn_V);

    if leng_V<MaxNumChanges_V
        turn_V(leng_V+1:MaxNumChanges_V)=turn_V(leng_V);
        warning('Too many turning points.')
    end

    turn_V_four(1:length(turn_V),i)=use_S_TP(turn_V);

end

%%
font_Size=14;F_Position=[2100,250,700,500];

f1a=figure

if S_bou~=0
    disp(['You choose the boundary restriction for S: ' num2str(S_bou)])
    he_gap1=area([0,Max_S],[S_bou 1-2*S_bou S_bou;S_bou 1-2*S_bou S_bou],0);
    he_gap1(2).FaceColor = [1,1,1];
    he_gap1(2).EdgeColor = [1,1,1];

    he_gap1(1).FaceColor = [0.4,0.4,0.4];
    he_gap1(1).EdgeColor = [0.4,0.4,0.4];
    he_gap1(1).FaceAlpha = 0.4;
    he_gap1(1).EdgeAlpha = 0.7;
    he_gap1(3).FaceColor = [0.4,0.4,0.4];
    he_gap1(3).EdgeColor = [0.4,0.4,0.4];
    he_gap1(3).FaceAlpha = 0.4;
    he_gap1(3).EdgeAlpha = 0.7;

    hold on
end

p1=plot(1:Max_S,V(1:Max_S),'b','LineWidth',2);
hold on

p2=plot([0 turn_V_four(1,used_sta_V)],[mean(V(1:turn_V_four(1,used_sta_V))) mean(V(1:turn_V_four(1,used_sta_V)))],...
    'Color','k','LineWidth',1.5,'LineStyle','--');
hold on
for i=1:MaxNumChanges_V-1
    plot([turn_V_four(i,used_sta_V) turn_V_four(i,used_sta_V)],[0 1],...
        'Color',[0.8 0.8 0.8],'LineWidth',1)
    hold on
    plot([turn_V_four(i,used_sta_V) turn_V_four(i+1,used_sta_V)],[mean(V(turn_V_four(i,used_sta_V):turn_V_four(i+1,used_sta_V))) mean(V(turn_V_four(i,used_sta_V):turn_V_four(i+1,used_sta_V)))],...
        'Color','k','LineWidth',1.5,'LineStyle','--');
    hold on
end
plot([turn_V_four(MaxNumChanges_V,used_sta_V) turn_V_four(MaxNumChanges_V,used_sta_V)],[0 1],...
    'Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
plot([turn_V_four(MaxNumChanges_V,used_sta_V) length(V)],[mean(V(turn_V_four(MaxNumChanges_V,used_sta_V):length(V))) mean(V(turn_V_four(MaxNumChanges_V,used_sta_V):length(V)))],...
    'Color','k','LineWidth',1.5,'LineStyle','--');
hold on

p3=plot(turn_V_four(:,used_sta_V),V(turn_V_four(:,used_sta_V)),'o','MarkerSize',5,...
    'MarkerFaceColor','r','MarkerEdgeColor','k');
hold on

xlim([0,Max_S])

ylabel('Concentration ratio (\lambda)','FontWeight','bold')
xlabel('Spherical Slepian functions (SSF)','FontWeight','bold')
legend([p3],{'TP(S)'}, ...
    'NumColumns',2,'Position',[0.5,0.84,0.4,0.03],'FontSize',font_Size+2,'box','off')

set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

% set(gca, 'xticklabel', get(gca, 'xtick'))
xlim([0,Max_S])

set(gcf,'Position',F_Position)

tif_name1a=[Attach_ALL 'TurningPoint_V.tif'];
print(f1a,'-dtiff','-r300',fullfile(fig_path_ALL,tif_name1a));

%% MSSA part

%% Should we fill the gap?

intit_num=size(A5_Centers_Sort,2)-1;

% for the each Slepian functions
[~,m]=size(A5_Centers_Sort(1).fill_deltacoffs);

%%%%%% The determination of (M, N) for gap-filling procedure is by other simulation

if any(A5_Centers_Sort(1).missing_months)

    disp(['Be careful: Here you need to choose the sliding window side M and cutoff number N for the gap-filling procedure.'])
    disp(['Be careful: The default set is: M = 96, N=8. Please be sure this is suitable for your length of data (M should not be larger than half of your data length)'])
    pause

    % consistent M
    M_gap=12*8; % 96 months

    % cutoff number N
    N_gap=8;

    disp(['The gap-filling procedure parameter: M = ' num2str(M_gap) ', N = ' num2str(N_gap)])

else

    disp('No gaps in the data.');
end

%% Gap filling (M-SSA) following (Gauer et al., 2022) and calculate the SSF
addpath('./MSSA')
% load(fullfile(getenv('IFILES'),['Results_' parrent_path],['Sen_' Attach_ALL 'c11.mat']));

intit_num=size(A5_Centers_Sort,2)-1;

% for the each Slepian functions
[~,m]=size(A5_Centers_Sort(1).fill_deltacoffs);

%
MSLAf_Centers_gapfilling=zeros(numel(fill_dates),intit_num,S_ol);
for ss=1:S_ol

    ss

    NaNflag=0;

    A5_Centers_Sort_mn=A5_Centers_Sort(1:intit_num);

    % substitute the MSSA_TS as the fill_deltacoffs for
    % each grid
    for ins=1:intit_num

        if ss<S+1
            A5_XXXX_mn=A5_Centers_Sort_mn(ins).fill_deltacoffs;
            A5_Centers_Sort_mn(ins).MSSA_TS=A5_XXXX_mn(:,ss)';
            A5_Centers_Sort_mn(ins).MSSA_TS_order=ss;
        else
            disp('we need to consider IB for ocean study.')
            A5_XXXX_mn=A5_Centers_Sort_mn(ins).fill_IBdeltacoffs;
            A5_Centers_Sort_mn(ins).MSSA_TS=A5_XXXX_mn';
            A5_Centers_Sort_mn(ins).MSSA_TS_order='IB';
        end

        if sum(isnan(A5_Centers_Sort_mn(ins).MSSA_TS))
            NaNflag=1;
            break
        end

    end

    % don't do MSSA for NaN
    if NaNflag
        warning(['NaN exists for m: ' num2str(ss)])
        continue
    end

    % replace the gap filling and do the final MSSA
    MSLAf_Centers=zeros(numel(fill_dates),intit_num);
    for ins=1:intit_num
        MSLAf_Centers(:,ins)=A5_Centers_Sort_mn(ins).MSSA_TS';
    end

    if any(A5_Centers_Sort(1).missing_months)
        % old version 1.0 (combine gap fitting and final, but may smooth some signals by the first replacement of 3 sigma)
        %         if mod(ss-1,5)==0
        %             [MSSAn_Total_CGJI_fill_reconst,MSSAn_RC,iter_num,chi,RCTest]=MSSA_gap_fitting(A5_Centers_Sort_mn,floor(numel(fill_dates)/2),Nc,...
        %                 fullfile(sub_fig_path_ALL,['MSSA_N' num2str(ss)]));
        %         else
        %             [MSSAn_Total_CGJI_fill_reconst,MSSAn_RC,iter_num,chi,RCTest]=MSSA_gap_fitting(A5_Centers_Sort_mn,floor(numel(fill_dates)/2),Nc);
        %
        %         end
        % MSSAf_Total_CGJI_fill_reconst=MSSAn_Total_CGJI_fill_reconst{max(iter_num)+1};

        % new version 1.1, only leakage values was used, without replacement of
        % 3 sigma
        %     if mod(ss-1,5)==0
        %         [MSSAn_CGJI_fill_reconst,~,iter_num,chi,~]=MSSA_gap_fitting_only(A5_Centers_Sort_mn,Mnum_rep,Nc1,...
        %             fullfile(fig_path_ALL,['MSSA_N' num2str(ss) '_Final']));
        %     else
        %         [MSSAn_CGJI_fill_reconst,~,iter_num,chi,~]=MSSA_gap_fitting_only(A5_Centers_Sort_mn,Mnum_rep,Nc1);
        %
        %     end
        %     replace_Centers=MSSAn_CGJI_fill_reconst{max(iter_num)+1};

        % new version 1.2, correct for calculating std only for missing values
        % 3 sigma
        if mod(ss-1,5)==0
            [MSSAn_CGJI_fill_reconst,~,iter_num,chi,~]=MSSA_gap_fitting_correct(A5_Centers_Sort_mn,M_gap,N_gap,...
                fullfile(fig_path_ALL,['MSSA_N' num2str(ss) '_Final']));
        else
            [MSSAn_CGJI_fill_reconst,~,iter_num,chi,~]=MSSA_gap_fitting_correct(A5_Centers_Sort_mn,M_gap,N_gap);

        end
        replace_Centers=MSSAn_CGJI_fill_reconst{max(iter_num)+1};

        for ins=1:intit_num

            MSLAf_Centers(A5_Centers_Sort(ins).missing_months,ins)=replace_Centers(A5_Centers_Sort(ins).missing_months,ins);

        end

    else

        iter_num=0;

    end

    MSLAf_Centers_gapfilling(:,:,ss)=MSLAf_Centers;

end

%% We need to first decide the optimal M for this region, and also calculate the turning point of N

%
disp('Now we want to search for the optimal M for reconstruction.');
disp('Notice that your data lengths had better exceeds 72, or the code will report error.')
disp('Or you could adjust the search M by adjust the variable ''test_M_rec''.')
pause
close all

A5_Centers_Sort_00=A5_Centers_Sort;
A5_Centers_Fit_00=A5_Centers_Fit;
MSLAf_Centers_gapfilling_00=MSLAf_Centers_gapfilling;
S_ol_00=S_ol;

%% determine the suitable M by the separation number
test_M_rec=[36:12:numel(fill_dates)/2];

RC_Seq=[];RC_RMS_Seq=[];
for i=1:length(test_M_rec)

    disp(['M = ' num2str(test_M_rec(i))]);

    for ss=1:S_ol

        for ins=1:intit_num
            if ss<S+1
                A5_Centers_Sort_mn(ins).MSSA_TS_order=ss;
            else
                disp('we need to consider IB for ocean study.')
                A5_Centers_Sort_mn(ins).MSSA_TS_order='IB';
            end
        end

        if intit_num*test_M_rec(i)<20
            warning('You used too small time series.')
            test_nssa=intit_num*test_M_rec(i);
        else
            test_nssa=20;
        end
        %             if mod(ss-1,5)==0
        %                 [MSSAf_CGJI_fill_reconst,MSSAn_RC,MSSAf_evalues,RCTest]=MSSA_final_noCDF(A5_Centers_Sort_mn,MSLAf_Centers_gapfilling(:,:,ss),Mssa2,N_rec,...
        %                     fullfile(fig_path_ALL,['MSSA_N' num2str(ss) '_Final']));
        %             else
        %     [MSSAf_CGJI_fill_reconst,MSSAn_RC,MSSAf_evalues,RCTest]=MSSA_final_noCDF(A5_Centers_Sort_mn,MSLAf_Centers_gapfilling(:,:,ss),Mssa2,N_rec(ss));
        [MSSAf_CGJI_fill_reconst,MSSAn_RC,MSSAf_evalues,RCTest]=MSSA_final_noCDF_freqsort(A5_Centers_Sort_mn,MSLAf_Centers_gapfilling(:,:,ss),test_M_rec(i),test_nssa);
        %             end

        % output
        %         CGJI_fill_MSSAcoffs(:,:,ss)=MSSAf_CGJI_fill_reconst;

        %         for j=1:size(MSSAf_evalues,1)
        %             MSSAf_evalues_sumup(j,:)=sum(MSSAf_evalues(1:j,:),1);
        %         end
        %         keep_mode=MSSAf_evalues_sumup(:,i)<Nssa_varupto;
        %         MSSAf_evalues_sumvariance(i)=sum(keep_mode);

        RC_Seq(i,ss,:)=mean(RCTest(5:7,:)');
        RC_RMS_Seq(i,ss,:)=mean(RCTest(5:7,:)'.*rms(mean(MSSAn_RC,3))');
    end

end


%%
Seq12_13_23=zeros(length(test_M_rec),3);
WSeq12_13_23=zeros(length(test_M_rec),3);

for i=1:length(test_M_rec)
    lg_str{i}=num2str(test_M_rec(i));
end

figure
subplot(2,3,1)
for i=1:length(test_M_rec)
    plot(1:S,squeeze(RC_Seq(i,1:S,1)));
    hold on
    Seq12_13_23(i,1)=mean(squeeze(RC_Seq(i,1:S,1)));
end
title('First Separation number (1/2)')
legend(lg_str)

subplot(2,3,2)
for i=1:length(test_M_rec)
    plot(1:S,squeeze(RC_Seq(i,1:S,2)));
    hold on
    Seq12_13_23(i,2)=mean(squeeze(RC_Seq(i,1:S,2)));
end
title('Second Separation number (1/3)')
legend(lg_str)

subplot(2,3,3)
for i=1:length(test_M_rec)
    plot(1:S,squeeze(RC_Seq(i,1:S,3)));
    hold on
    Seq12_13_23(i,3)=mean(squeeze(RC_Seq(i,1:S,3)));
end
title('Third Separation number (2/3)')
legend(lg_str)

subplot(2,3,4)
for i=1:length(test_M_rec)
    %     plot(1:N_ol,squeeze(RC_Seq(i,:,1))./N_rec);
    %     hold on
    %     WSeq12_13_23(i,1)=mean(squeeze(RC_Seq(i,:,1))./N_rec);

    plot(1:S,squeeze(RC_RMS_Seq(i,1:S,1)));
    hold on
    WSeq12_13_23(i,1)=mean(squeeze(RC_RMS_Seq(i,1:S,1)));

end
title('First Separation number (1/2)')
legend(lg_str)

subplot(2,3,5)
for i=1:length(test_M_rec)
    %     plot(1:N_ol,squeeze(RC_Seq(i,:,2))./N_rec);
    %     hold on
    %     WSeq12_13_23(i,2)=mean(squeeze(RC_Seq(i,:,2))./N_rec);

    plot(1:S,squeeze(RC_RMS_Seq(i,1:S,2)));
    hold on
    WSeq12_13_23(i,2)=mean(squeeze(RC_RMS_Seq(i,1:S,2)));
end
title('Second Separation number (1/3)')
legend(lg_str)

subplot(2,3,6)
for i=1:length(test_M_rec)
    %     plot(1:N_ol,squeeze(RC_Seq(i,:,3))./N_rec);
    %     hold on
    %     WSeq12_13_23(i,3)=mean(squeeze(RC_Seq(i,:,3))./N_rec);

    plot(1:S,squeeze(RC_RMS_Seq(i,1:S,3)));
    hold on
    WSeq12_13_23(i,3)=mean(squeeze(RC_RMS_Seq(i,1:S,3)));
end
title('Third Separation number (2/3)')
legend(lg_str)

% simple version, use the maximum to determine the M
[WSeq_max,WSeq_num]=max(WSeq12_13_23);
test_M_rec_use=WSeq_num(1);
M_rec=test_M_rec(WSeq_num(1));

% breaking point version, use the break point
% used_sta_M=4; % what kind of turning point do you want? (default is RSS)
% MaxNumChanges_M=1; % You only need one optimal M
%
% option_sta=["mean","rms","std","linear"];
% turn_M_four=[];turn_WM_four=[];
% for i=1:4
%     used_sta=char(option_sta(i));
%     % whole SCS
%     [turn_M] = findchangepts(Seq12_13_23(:,1), 'Statistic', used_sta, 'MaxNumChanges', MaxNumChanges_M);
%     [turn_WM] = findchangepts(WSeq12_13_23(:,1), 'Statistic', used_sta, 'MaxNumChanges', MaxNumChanges_M);
%
%     turn_M_four(:,i)=turn_M;
%     turn_WM_four(:,i)=turn_WM;
% end

% test_M_rec_use=turn_WM_four(:,used_sta_M);
% M_rec=test_M_rec(turn_WM_four(:,used_sta_M));

%% plot for the optimal M
M_color=mymap('tab10');
font_Size=14;F_Position=[2100,250,1000,450];
F_gap=[.03 .08];F_marg_h=[.13 .05];F_marg_w=[.07 .02];

fm=figure

ha=tight_subplot(1,2,F_gap,F_marg_h,F_marg_w);

set(gcf,'Position',F_Position)

axes(ha(1))
for i=1:length(test_M_rec)
    %     plot(1:N_ol,squeeze(RC_Seq(i,:,1))./N_rec);
    %     hold on
    %     WSeq12_13_23(i,1)=mean(squeeze(RC_Seq(i,:,1))./N_rec);

    plot(1:S,squeeze(RC_RMS_Seq(i,1:S,1)),'color',M_color(26*i-25,:),...
        'LineWidth',1.2);
    hold on
    WSeq12_13_23(i,1)=mean(squeeze(RC_RMS_Seq(i,1:S,1)));

end

xlim([0,S+1])
legend(lg_str,'box','off','NumColumns',2,'FontSize',font_Size-2)
ylabel('Separation index','FontWeight','bold')
xlabel('Spherical Slepian coefficients (SSF)','FontWeight','bold')

% title(Area)

set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

axes(ha(2))
plot(test_M_rec,WSeq12_13_23(:,1),'LineWidth',2,'color','k');
hold on
% p1=plot(M_rec,WSeq12_13_23(turn_WM_four(:,used_sta_M),1),'o','MarkerSize',8,...
%     'MarkerFaceColor','r','MarkerEdgeColor','k');
p1=plot(M_rec,WSeq12_13_23(test_M_rec_use,1),'o','MarkerSize',8,...
    'MarkerFaceColor','r','MarkerEdgeColor','k');

xticks(test_M_rec)

ylabel('Mean separation index of the first 50 SSF','FontWeight','bold')
xlabel('Size of Sliding window M (month)','FontWeight','bold')

legend(p1,'Optimal M','box','off','Location','northwest','FontSize',font_Size)
xlim([test_M_rec(1)-5,test_M_rec(end)+5])
set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

tif_namefm=[Attach_ALL 'M_rec_Determine.tif'];
print(fm,'-dtiff','-r300',fullfile(fig_path_ALL,tif_namefm));

disp(['So the optimal M for reconstruction is determined to be : ' num2str(M_rec) ])

save('MSSA_Information.mat','M_rec','-append')

pause(0.5)
%% Slepian

if strcmp(A5_Centers_Sort_00(1).process{3},'ocean')
    Max_S_ol=Max_S+1;
else
    Max_S_ol=Max_S;
end

% for MSSA

MSSAf_evalues_mm=zeros(intit_num*M_rec,Max_S_ol);
for ss=1:Max_S_ol
    [~,~,MSSAf_evalues,~]=MSSA_final_noCDF_freqsort(A5_Centers_Sort_mn,MSLAf_Centers_gapfilling(:,:,ss),M_rec,N_gap); %error Mssa1 should be Mssa2 to be consistent with the calculation of M*4

    MSSAf_evalues_mm(:,ss)=MSSAf_evalues;
end

MSSAf_evalues_sumup=zeros(size(MSSAf_evalues_mm));
for i=1:size(MSSAf_evalues_mm,1)
    MSSAf_evalues_sumup(i,:)=sum(MSSAf_evalues_mm(1:i,:),1);
end


%% turning point of the sum of eigenvalues of MSSA (MSSA part)
option_sta=["mean","rms","std","linear"];

used_sta_MSSA=4; % what kind of turning point do you want? (default is RSS)
MaxNumChanges_MSSA=Turning_number; % how many turning point do you want?
% figure
% subplot(211)
% plot(1:size(MSSAf_evalues_sumup,2),MSSAf_evalues_sumup)
% subplot(212)
% plot(1:size(MSSAf_evalues_sumup,1),MSSAf_evalues_sumup')

turn_MSSA_four={};
for j=1:size(MSSAf_evalues_sumup,2)
    turn_MSSA_four{j}=nan(MaxNumChanges_MSSA,4);
    for i=1:4
        used_sta=char(option_sta(i));
        %         turn_MSSA_four{j}=nan(max_turning,4);

        use_N_TP=find(MSSAf_evalues_sumup(:,j)>threshold_N(1) & MSSAf_evalues_sumup(:,j)<threshold_N(2));

        % whole SCS
        [turn_MSSA] = findchangepts(MSSAf_evalues_sumup(use_N_TP,j), 'Statistic', used_sta, 'MaxNumChanges', MaxNumChanges_MSSA);

        leng_MSSA=length(turn_MSSA);

        if leng_MSSA<MaxNumChanges_MSSA
            turn_MSSA(leng_MSSA+1:MaxNumChanges_MSSA)=turn_MSSA(leng_MSSA);
            warning('Too many turning points for MSSA.')
        end

        turn_MSSA_four{j}(1:length(turn_MSSA),i)=use_N_TP(turn_MSSA);


    end
end

% MSSAf_evalues_sumvariance=zeros(1,N_ol);
% for i=1:N_ol
%     keep_mode=MSSAf_evalues_sumup(:,i)<Nssa_varupto;
%     MSSAf_evalues_sumvariance(i)=sum(keep_mode);
%
% %     Nssa2_SN(:,i)=turn_MSSA_four{j}
% end

% N_rec
% N_range(nn)=round(mean(MSSAf_evalues_sumvariance(1:S))); %do not consider the IB here
%
% N_rec=N_range(nn);
%
% if N_rec==0
%     N_rec=1; % at least remain one RC
% end

%% turning point MSSA
font_Size=14;F_Position=[2100,250,700,500];

f1b=figure

max_mode=min(150,intit_num*M_rec);

if N_bou~=0
    disp(['You choose the boundary restriction for N: ' num2str(N_bou)])
    he_gap1=area([0,max_mode],[N_bou 1-2*N_bou N_bou;N_bou 1-2*N_bou N_bou],0);
    he_gap1(2).FaceColor = [1,1,1];
    he_gap1(2).EdgeColor = [1,1,1];

    he_gap1(1).FaceColor = [0.4,0.4,0.4];
    he_gap1(1).EdgeColor = [0.4,0.4,0.4];
    he_gap1(1).FaceAlpha = 0.4;
    he_gap1(1).EdgeAlpha = 0.7;
    he_gap1(3).FaceColor = [0.4,0.4,0.4];
    he_gap1(3).EdgeColor = [0.4,0.4,0.4];
    he_gap1(3).FaceAlpha = 0.4;
    he_gap1(3).EdgeAlpha = 0.7;

    hold on
end

colorr=mymap('Blues');
dd=linspace(1,256,Max_S);
map_colori=colorr(floor(flipud(dd')),:);TickLabels_c={};
for j=1:Max_S
    %     color_line=repmat(0.2+0.02*j,1,3);
    color_line=map_colori(j,:);
    plot(1:max_mode,MSSAf_evalues_sumup(1:max_mode,j),'color',color_line,'LineWidth',1);
    hold on

    p3=plot(turn_MSSA_four{j}(:,used_sta_MSSA),MSSAf_evalues_sumup(turn_MSSA_four{j}(:,used_sta_MSSA),j),'o','MarkerSize',5,...
        'MarkerFaceColor','y','MarkerEdgeColor','k');
    hold on

    xlim([0,max_mode])

    for i=1:MaxNumChanges_MSSA
        turning_MSSA(j,i)=turn_MSSA_four{j}(i,used_sta_MSSA);
        turning_MSSA_evalues(j,i)=MSSAf_evalues_sumup(turn_MSSA_four{j}(i,used_sta_MSSA),j);
    end

    TickLabels_c{j}=num2str(j);
end

colormap(gca,map_colori)
%         colorbar('Ticks',linspace(0,1,iter_num(i)),...
%             'TickLabels',TickLabels_c)
cb=colorbar('Ticks',[1/Max_S/2:10/Max_S:1,1],...
    'TickLabels',TickLabels_c([1:10:end,end]))

for i=1:MaxNumChanges_MSSA
    plot(turning_MSSA(1:Max_S,:),turning_MSSA_evalues(1:Max_S,:),'k')
end

% cb = colorbar;
cb.Position = [0.8 0.2 0.04,0.3];
title(cb,'SSF','fontsize',14);

ylim([0,1])

ylabel({'\bf Sum of normalized eigenvalues ($\sum\widetilde{\lambda}$)','\bf of M-SSA'},...
    'Interpreter','latex','FontWeight','bold')
xlabel('Reconstructed components (RC) of each SSF','FontWeight','bold')

legend([p3],{'TP(N)'}, ...
    'NumColumns',2,'Position',[0.55,0.7,0.4,0.03],'FontSize',font_Size+2,'box','off')

set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

xlim([0,max_mode])

% set(gca, 'xticklabel', get(gca, 'xtick'))
% xlim([0,Max_SSF])

set(gcf,'Position',F_Position)

tif_name1b=[Attach_ALL 'TurningPoint_MSSA.tif'];
print(f1b,'-dtiff','-r300',fullfile(fig_path_ALL,tif_name1b));

%% combination of V and MSSA
font_Size=14;F_Position=[2100,250,1400,500];
F_gap=[.03 .08];F_marg_h=[.12 .05];F_marg_w=[.07 .02];

f1c=figure

set(gcf,'Position',F_Position)

ha=tight_subplot(1,2,F_gap,F_marg_h,F_marg_w);

axes(ha(1))

if S_bou~=0
    disp(['You choose the boundary restriction for S: ' num2str(S_bou)])
    he_gap1=area([0,Max_S],[S_bou 1-2*S_bou S_bou;S_bou 1-2*S_bou S_bou],0);
    he_gap1(2).FaceColor = [1,1,1];
    he_gap1(2).EdgeColor = [1,1,1];

    he_gap1(1).FaceColor = [0.4,0.4,0.4];
    he_gap1(1).EdgeColor = [0.4,0.4,0.4];
    he_gap1(1).FaceAlpha = 0.4;
    he_gap1(1).EdgeAlpha = 0.7;
    he_gap1(3).FaceColor = [0.4,0.4,0.4];
    he_gap1(3).EdgeColor = [0.4,0.4,0.4];
    he_gap1(3).FaceAlpha = 0.4;
    he_gap1(3).EdgeAlpha = 0.7;

    hold on
end

p1=plot(1:Max_S,V(1:Max_S),'Color',[0 0.4470 0.7410],'LineWidth',2);
hold on
for i=1:MaxNumChanges_V
    text(turn_V_four(i,used_sta_V)+2,V(turn_V_four(i,used_sta_V)),['S(' num2str(i) ')'],...
        'FontSize',font_Size,'Color','r')
    hold on
end

p2=plot([0 turn_V_four(1,used_sta_V)],[mean(V(1:turn_V_four(1,used_sta_V))) mean(V(1:turn_V_four(1,used_sta_V)))],...
    'Color','k','LineWidth',1.5,'LineStyle','--');
hold on
for i=1:MaxNumChanges_V-1
    plot([turn_V_four(i,used_sta_V) turn_V_four(i,used_sta_V)],[0 1],...
        'Color',[0.8 0.8 0.8],'LineWidth',1)
    hold on
    plot([turn_V_four(i,used_sta_V) turn_V_four(i+1,used_sta_V)],[mean(V(turn_V_four(i,used_sta_V):turn_V_four(i+1,used_sta_V))) mean(V(turn_V_four(i,used_sta_V):turn_V_four(i+1,used_sta_V)))],...
        'Color','k','LineWidth',1.5,'LineStyle','--');
    hold on
end
plot([turn_V_four(MaxNumChanges_V,used_sta_V) turn_V_four(MaxNumChanges_V,used_sta_V)],[0 1],...
    'Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
plot([turn_V_four(MaxNumChanges_V,used_sta_V) length(V)],[mean(V(turn_V_four(MaxNumChanges_V,used_sta_V):length(V))) mean(V(turn_V_four(MaxNumChanges_V,used_sta_V):length(V)))],...
    'Color','k','LineWidth',1.5,'LineStyle','--');
hold on

p3=plot(turn_V_four(:,used_sta_V),V(turn_V_four(:,used_sta_V)),'o','MarkerSize',5,...
    'MarkerFaceColor','r','MarkerEdgeColor','k');
hold on

xlim([0,Max_S])

ylabel({'Concentration ratio (\lambda)','of Slepian conversion'},'FontWeight','bold')
xlabel('Spherical Slepian Coefficients (SSC)','FontWeight','bold')
legend([p3],{['TP-S(1~' num2str(MaxNumChanges_V) ')']}, ...
    'NumColumns',2,'Position',[0.3,0.7,0.2,0.1],'FontSize',font_Size+2,'box','off')

set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

% set(gca, 'xticklabel', get(gca, 'xtick'))
xlim([0,Max_S])

axes(ha(2))
% f1b=figure

if N_bou~=0
    disp(['You choose the boundary restriction for N: ' num2str(N_bou)])
    he_gap1=area([0,max_mode],[N_bou 1-2*N_bou N_bou;N_bou 1-2*N_bou N_bou],0);
    he_gap1(2).FaceColor = [1,1,1];
    he_gap1(2).EdgeColor = [1,1,1];

    he_gap1(1).FaceColor = [0.4,0.4,0.4];
    he_gap1(1).EdgeColor = [0.4,0.4,0.4];
    he_gap1(1).FaceAlpha = 0.4;
    he_gap1(1).EdgeAlpha = 0.7;
    he_gap1(3).FaceColor = [0.4,0.4,0.4];
    he_gap1(3).EdgeColor = [0.4,0.4,0.4];
    he_gap1(3).FaceAlpha = 0.4;
    he_gap1(3).EdgeAlpha = 0.7;

    hold on
end

colorr=mymap('Blues');
dd=linspace(1,256,Max_S);
map_colori=colorr(floor(flipud(dd')),:);
max_mode=min(150,intit_num*M_rec);TickLabels_c={};
for j=1:Max_S
    %     color_line=repmat(0.2+0.02*j,1,3);
    color_line=map_colori(j,:);
    plot(1:max_mode,MSSAf_evalues_sumup(1:max_mode,j),'color',color_line,'LineWidth',1);
    hold on

    p3=plot(turn_MSSA_four{j}(:,used_sta_MSSA),MSSAf_evalues_sumup(turn_MSSA_four{j}(:,used_sta_MSSA),j),'o','MarkerSize',5,...
        'MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerEdgeColor','k');
    hold on

    xlim([0,max_mode])

    for i=1:MaxNumChanges_MSSA
        turning_MSSA(j,i)=turn_MSSA_four{j}(i,used_sta_MSSA);
        turning_MSSA_evalues(j,i)=MSSAf_evalues_sumup(turn_MSSA_four{j}(i,used_sta_MSSA),j);
    end

    TickLabels_c{j}=num2str(j);
end

colormap(gca,map_colori)
%         colorbar('Ticks',linspace(0,1,iter_num(i)),...
%             'TickLabels',TickLabels_c)
cb=colorbar('Ticks',[1/Max_S/2:10/Max_S:1,1],...
    'TickLabels',TickLabels_c([1:10:end,end]))

for i=1:MaxNumChanges_MSSA
    plot(turning_MSSA(1:Max_S,:),turning_MSSA_evalues(1:Max_S,:),'k')

    hold on
    text(max(turning_MSSA(:,i))+1,min(turning_MSSA_evalues(:,i))-0.01,['N(' num2str(i) ')'],...
        'FontSize',font_Size,'Color',[0.9290 0.6940 0.1250])
    hold on
end

% cb = colorbar;
cb.Position = [0.9 0.2 0.02,0.3];
title(cb,'SSF','fontsize',14);

% xlim([0,max_mode])

ylabel({'\bf Sum of normalized eigenvalues ($\sum\widetilde{\lambda}$)','\bf of M-SSA'},...
    'Interpreter','latex','FontWeight','bold')
xlabel('Reconstructed components (RC) of each SSC','FontWeight','bold')

legend([p3],{['TP-N(1~' num2str(MaxNumChanges_MSSA) ')']}, ...
    'NumColumns',2,'Position',[0.8,0.7,0.2,0.1],'FontSize',font_Size+2,'box','off')

set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

xlim([0,max_mode])
ylim([0,1])
% set(gca, 'xticklabel', get(gca, 'xtick'))
% xlim([0,Max_SSF])

% set(gcf,'Position',F_Position)

tif_name1c=['TurningPoint_V_MSSA.tif'];
print(f1c,'-dtiff','-r300',fullfile(fig_path_ALL,tif_name1c));

save(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' Attach_ALL ...
    '_TP' num2str(Turning_number) '_M' num2str(M_rec) '_Breakingpoint.mat']),...
    'MaxNumChanges_MSSA','MaxNumChanges_V','turn_V_four','turn_MSSA_four',...
    'N_gap','M_gap','M_rec','Turning_number','RC_RMS_Seq','V',...
    'used_sta_V','used_sta_MSSA','MSSAf_evalues_sumup','Max_S','test_M_rec',...
    'threshold_S','threshold_N')

% figure
% plot(XY(:,1),XY(:,2))
% hold on
% plot(XY_ori(:,1),XY_ori(:,2))

disp('Now we have already determined the breaking points of S and N. ')
disp('Then, we need the remove the noise by statistical method at different significance level p.')
pause
%% for different S and N for reconstruction
SIG_use=[0.05,0.1,0.3];

for Noise_SigLev=SIG_use

    % the TP(S1~Sn)
    for TP_ss=1:MaxNumChanges_V

        %     for TP_ss=MaxNumChanges_V
        %       for  TP_ss=5;
        %           for  TP_ss=S_shannon

        % the TP(N1~Nn)
        for TP_nn=1:MaxNumChanges_MSSA

            %        for  TP_nn=MaxNumChanges_MSSA
            %             for    TP_nn=5

            % initialize
            A5_Centers_Sort=A5_Centers_Sort_00;
            A5_Centers_Fit=A5_Centers_Fit_00;

            % for slepian
            S=turn_V_four(TP_ss,used_sta_V);
            % only for shannon
            %         S=TP_ss;
            S_sig=S;

            N_rec=[];
            for k=1:S
                N_rec(k)=turn_MSSA_four{k}(TP_nn,used_sta_MSSA);
            end

            if strcmp(A5_Centers_Sort(1).process{3},'ocean')
                S_ol=S+1;
                N_rec(S+1)=turn_MSSA_four{S_ol_00}(TP_nn,used_sta_MSSA);
            else
                S_ol=S;
            end

            MSLAf_Centers_gapfilling=MSLAf_Centers_gapfilling_00;
            MSLAf_Centers_gapfilling(:,:,S_ol)=MSLAf_Centers_gapfilling_00(:,:,S_ol_00);

            %         N_rec=max_mode;

            %% reconstruction (M-SSA) following (Gauer et al., 2022)

            CGJI_fill_MSSAcoffs=nan(numel(fill_dates),intit_num,S_ol);
            CGJI_fill_MSSAcoffs_iter_num=nan(intit_num,S_ol);
            CGJI_fill_MSSAcoffs_evalues=[];
            CGJI_fill_MSSAcoffs_RC=cell(S_ol,1);
            CGJI_fill_MSSAcoffs_RCTest=cell(S_ol,1);

            for ss=1:S_ol

                for ins=1:intit_num
                    if ss<S+1
                        A5_Centers_Sort_mn(ins).MSSA_TS_order=ss;
                    else
                        disp('we need to consider IB for ocean study.')
                        A5_Centers_Sort_mn(ins).MSSA_TS_order='IB';
                    end
                end

                %             if mod(ss-1,5)==0
                %                 [MSSAf_CGJI_fill_reconst,MSSAn_RC,MSSAf_evalues,RCTest]=MSSA_final_noCDF(A5_Centers_Sort_mn,MSLAf_Centers_gapfilling(:,:,ss),M_rec,N_rec,...
                %                     fullfile(fig_path_ALL,['MSSA_N' num2str(ss) '_Final']));
                %             else
                %             [MSSAf_CGJI_fill_reconst,MSSAn_RC,MSSAf_evalues,RCTest]=MSSA_final_noCDF(A5_Centers_Sort_mn,MSLAf_Centers_gapfilling(:,:,ss),M_rec,N_rec(ss));
                [MSSAf_CGJI_fill_reconst,MSSAn_RC,MSSAf_evalues,RCTest]=MSSA_final_noCDF_freqsort(A5_Centers_Sort_mn,MSLAf_Centers_gapfilling(:,:,ss),M_rec,N_rec(ss),...
                    [],Noise_SigLev);

                %             end

                % output
                %         CGJI_fill_MSSAcoffs(:,:,ss)=MSSAf_CGJI_fill_reconst;
                CGJI_fill_MSSAcoffs(:,:,ss)=MSSAf_CGJI_fill_reconst;
                CGJI_fill_MSSAcoffs_iter_num(:,ss)=iter_num;
                CGJI_fill_MSSAcoffs_evalues(:,ss)=MSSAf_evalues;
                CGJI_fill_MSSAcoffs_RC{ss}=MSSAn_RC;
                CGJI_fill_MSSAcoffs_RCTest{ss}=RCTest;
            end

            %%
            % calculate the mean values
            % criteria for valid data
            % reject the lillietest test and K-S test to meet the CDF test
            CGJI_fill_MSSAcoffs_RC_both=cell(S_ol,1);
            CGJI_fill_MSSAcoffs_RC_noise=cell(S_ol,1);
            CGJI_fill_MSSAcoffs_eva=[];

            RCFreDomi=cell(S_ol,1);
            RC_noise=cell(S_ol,1);
            RC_Seq=zeros(S_ol,3);
            RC_p=cell(S_ol,1);

            CGJI_fill_MSSAcoffs_both=zeros(numel(fill_dates),intit_num,S_ol);
            CGJI_fill_MSSAcoffs_noise=zeros(numel(fill_dates),intit_num,S_ol);
            CGJI_fill_MSSAcoffs_resid=zeros(numel(fill_dates),intit_num,S_ol);

            CGJI_fill_MSSAcoffs_both_annual=zeros(numel(fill_dates),intit_num,S_ol);
            CGJI_fill_MSSAcoffs_both_long=zeros(numel(fill_dates),intit_num,S_ol);
            CGJI_fill_MSSAcoffs_both_semi=zeros(numel(fill_dates),intit_num,S_ol);
            CGJI_fill_MSSAcoffs_both_short=zeros(numel(fill_dates),intit_num,S_ol);

            CGJI_fill_MSSAcoffs_noise_annual=zeros(numel(fill_dates),intit_num,S_ol);
            CGJI_fill_MSSAcoffs_noise_long=zeros(numel(fill_dates),intit_num,S_ol);
            CGJI_fill_MSSAcoffs_noise_semi=zeros(numel(fill_dates),intit_num,S_ol);
            CGJI_fill_MSSAcoffs_noise_short=zeros(numel(fill_dates),intit_num,S_ol);

            for ss=1:S_ol
                RCTest= CGJI_fill_MSSAcoffs_RCTest{ss};
                valid_both=find(RCTest(1,:)+RCTest(2,:)==2);
                valid_noise=find(RCTest(1,:)+RCTest(2,:)<2);

                RCFreDomi {ss}=  CGJI_fill_MSSAcoffs_RCTest{ss}(4,:)';
                %             RC_p {ss}=  RCTest(8:9,:)';
                %             RC_noise {ss} = zeros(N_rec(ss),1);
                %             RC_noise {ss}(valid_noise)=1;

                valid_both_annual=find(RCTest(1,:)+RCTest(2,:)==2 & RCTest(4,:)==1);
                valid_both_long=find(RCTest(1,:)+RCTest(2,:)==2 & RCTest(4,:)<1);
                valid_both_semi=find(RCTest(1,:)+RCTest(2,:)==2 & RCTest(4,:)==2);
                %     valid_both_short=setdiff(valid_both,[valid_both_annual,valid_both_long,valid_both_semi]); % notice that short period should include the semi-annual
                valid_both_short=setdiff(valid_both,[valid_both_annual,valid_both_long]); % this is the corrected one

                valid_noise_annual=find(RCTest(1,:)+RCTest(2,:)<2 & RCTest(4,:)==1);
                valid_noise_long=find(RCTest(1,:)+RCTest(2,:)<2 & RCTest(4,:)<1);
                valid_noise_semi=find(RCTest(1,:)+RCTest(2,:)<2 & RCTest(4,:)==2);
                %     valid_noise_short=setdiff(valid_noise,[valid_noise_annual,valid_noise_long,valid_noise_semi]); % notice that short period should include the semi-annual
                valid_noise_short=setdiff(valid_noise,[valid_noise_annual,valid_noise_long]); % this is the corrected one

                RC_Seq(ss,:)=mean(RCTest(5:7,:)');

                %     valid_both=valid_both(valid_both<=Nc2_threshold);

                CGJI_fill_MSSAcoffs_RC_both{ss,1}=valid_both;
                CGJI_fill_MSSAcoffs_RC_both{ss,2}=valid_both_long;
                CGJI_fill_MSSAcoffs_RC_both{ss,3}=valid_both_annual;
                CGJI_fill_MSSAcoffs_RC_both{ss,4}=valid_both_semi;
                CGJI_fill_MSSAcoffs_RC_both{ss,5}=valid_both_short;
                CGJI_fill_MSSAcoffs_RC_noise{ss,1}=valid_noise;
                CGJI_fill_MSSAcoffs_RC_noise{ss,2}=valid_noise_long;
                CGJI_fill_MSSAcoffs_RC_noise{ss,3}=valid_noise_annual;
                CGJI_fill_MSSAcoffs_RC_noise{ss,4}=valid_noise_semi;
                CGJI_fill_MSSAcoffs_RC_noise{ss,5}=valid_noise_short;

                MSSAn_RC=CGJI_fill_MSSAcoffs_RC{ss};

                CGJI_fill_MSSAcoffs_both(:,:,ss)=squeeze(sum(MSSAn_RC(:,valid_both,:),2));
                CGJI_fill_MSSAcoffs_both_annual(:,:,ss)=squeeze(sum(MSSAn_RC(:,valid_both_annual,:),2));
                CGJI_fill_MSSAcoffs_both_long(:,:,ss)=squeeze(sum(MSSAn_RC(:,valid_both_long,:),2));
                CGJI_fill_MSSAcoffs_both_semi(:,:,ss)=squeeze(sum(MSSAn_RC(:,valid_both_semi,:),2));
                CGJI_fill_MSSAcoffs_both_short(:,:,ss)=squeeze(sum(MSSAn_RC(:,valid_both_short,:),2));

                CGJI_fill_MSSAcoffs_noise(:,:,ss)=squeeze(sum(MSSAn_RC(:,valid_noise,:),2));
                %     CGJI_fill_MSSAcoffs_noise_ori(:,:,ss)=squeeze(sum(MSSAn_RC,2)-sum(MSSAn_RC(:,valid_both,:),2));     % this should be the same as the noise
                CGJI_fill_MSSAcoffs_noise_annual(:,:,ss)=squeeze(sum(MSSAn_RC(:,valid_noise_annual,:),2));
                CGJI_fill_MSSAcoffs_noise_long(:,:,ss)=squeeze(sum(MSSAn_RC(:,valid_noise_long,:),2));
                CGJI_fill_MSSAcoffs_noise_semi(:,:,ss)=squeeze(sum(MSSAn_RC(:,valid_noise_semi,:),2));
                CGJI_fill_MSSAcoffs_noise_short(:,:,ss)=squeeze(sum(MSSAn_RC(:,valid_noise_short,:),2));

                CGJI_fill_MSSAcoffs_resid(:,:,ss)=MSLAf_Centers_gapfilling(:,:,ss)-squeeze(sum(MSSAn_RC,2));
                CGJI_fill_MSSAcoffs_noise(:,:,ss)=squeeze(sum(MSSAn_RC,2)-sum(MSSAn_RC(:,valid_both,:),2));


                CGJI_fill_MSSAcoffs_eva(ss,:)=[sum(CGJI_fill_MSSAcoffs_evalues(1:N_rec(ss),ss),1),...
                    sum(CGJI_fill_MSSAcoffs_evalues(valid_both,ss),1),...
                    sum(CGJI_fill_MSSAcoffs_evalues(valid_both_long,ss),1),...
                    sum(CGJI_fill_MSSAcoffs_evalues(valid_both_annual,ss),1),...
                    sum(CGJI_fill_MSSAcoffs_evalues(valid_both_semi,ss),1),...
                    sum(CGJI_fill_MSSAcoffs_evalues(valid_both_short,ss),1),...
                    sum(CGJI_fill_MSSAcoffs_evalues(valid_noise,ss),1),...
                    sum(CGJI_fill_MSSAcoffs_evalues(valid_noise_long,ss),1),...
                    sum(CGJI_fill_MSSAcoffs_evalues(valid_noise_annual,ss),1),...
                    sum(CGJI_fill_MSSAcoffs_evalues(valid_noise_semi,ss),1),...
                    sum(CGJI_fill_MSSAcoffs_evalues(valid_noise_short,ss),1),...
                    ];
            end

            %%
            for ins=1:intit_num
                A5_Centers_Sort(ins).fill_MSSAcoffs_original=squeeze(MSLAf_Centers_gapfilling(:,ins,:));
                A5_Centers_Sort(ins).fill_MSSAcoffs=squeeze(CGJI_fill_MSSAcoffs(:,ins,:));
                A5_Centers_Sort(ins).fill_MSSAcoffs_IterationN=squeeze(CGJI_fill_MSSAcoffs_iter_num(ins,:));

                A5_Centers_Sort(ins).fill_MSSAcoffs_RC=squeeze(CGJI_fill_MSSAcoffs_RC);
                A5_Centers_Sort(ins).fill_MSSAcoffs_both=squeeze(CGJI_fill_MSSAcoffs_both(:,ins,:));
                A5_Centers_Sort(ins).fill_MSSAcoffs_both_annual=squeeze(CGJI_fill_MSSAcoffs_both_annual(:,ins,:));
                A5_Centers_Sort(ins).fill_MSSAcoffs_both_long=squeeze(CGJI_fill_MSSAcoffs_both_long(:,ins,:));
                A5_Centers_Sort(ins).fill_MSSAcoffs_both_semi=squeeze(CGJI_fill_MSSAcoffs_both_semi(:,ins,:));
                A5_Centers_Sort(ins).fill_MSSAcoffs_both_short=squeeze(CGJI_fill_MSSAcoffs_both_short(:,ins,:));

                A5_Centers_Sort(ins).fill_MSSAcoffs_noise=squeeze(CGJI_fill_MSSAcoffs_noise(:,ins,:));
                A5_Centers_Sort(ins).fill_MSSAcoffs_resid=squeeze(CGJI_fill_MSSAcoffs_resid(:,ins,:));

                A5_Centers_Sort(ins).fill_MSSAcoffs_eva=CGJI_fill_MSSAcoffs_eva;

                A5_Centers_Sort(ins).fill_MSSAcoffs_RC=cell(S_ol,1);
                for ss=1:S_ol
                    A5_Centers_Sort(ins).fill_MSSAcoffs_RC{ss}=squeeze(CGJI_fill_MSSAcoffs_RC{ss}(:,:,ins));
                end
            end

            A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_original=squeeze(sum(MSLAf_Centers_gapfilling,2)/intit_num);
            A5_Centers_Sort(intit_num+1).fill_MSSAcoffs=squeeze(sum(CGJI_fill_MSSAcoffs,2)/intit_num);
            A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_IterationN=squeeze(sum(CGJI_fill_MSSAcoffs_iter_num,1)/intit_num);

            % A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_signalRC=CGJI_fill_MSSA_RC_signal;
            A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_both=squeeze(sum(CGJI_fill_MSSAcoffs_both,2)/intit_num);
            A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_both_annual=squeeze(sum(CGJI_fill_MSSAcoffs_both_annual,2)/intit_num);
            A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_both_long=squeeze(sum(CGJI_fill_MSSAcoffs_both_long,2)/intit_num);
            A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_both_semi=squeeze(sum(CGJI_fill_MSSAcoffs_both_semi,2)/intit_num);
            A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_both_short=squeeze(sum(CGJI_fill_MSSAcoffs_both_short,2)/intit_num);

            A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_noise=squeeze(sum(CGJI_fill_MSSAcoffs_noise,2)/intit_num);
            A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_resid=squeeze(sum(CGJI_fill_MSSAcoffs_resid,2)/intit_num);

            A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_eva=CGJI_fill_MSSAcoffs_eva;

            A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_RC=cell(S_ol,1);
            for ss=1:S_ol
                RC_sum=zeros(size(A5_Centers_Sort(1).fill_MSSAcoffs_RC{ss}));
                for ins=1:intit_num
                    RC_sum=RC_sum+A5_Centers_Sort(ins).fill_MSSAcoffs_RC{ss};
                end
                A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_RC{ss}=RC_sum/intit_num;
            end

            A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_RCTest=CGJI_fill_MSSAcoffs_RCTest;

            close all

            %         save(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' Attach_ALL '_pro.mat']),...
            %             'A5_Centers_Sort', 'A5_Centers_Fit','c11cmn','Lwindow','S','S_sig','N_ol','S_shannon','CC','SHC_degord','BasinArea')

            %% calculate the explained eigenvalues for each coefficient for each components (e.g., annual)
            %         load(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' Attach_ALL '_pro.mat']));

            % if the longitude ranges from 0 to 360. Make it consistent to [-180,180]
            if XY(1,1)>180
                XY(:,1)=XY(:,1)-360;
            end

            if strcmp(parrent_path(5:end),'GL')
                c11cmn_comTo12_lat=1:length(lat);
                c11cmn_comTo12_lon=1:length(lon);
            else
                c11cmn_comTo12_lat=(c11cmn_combine(2)-c11cmn(2)+1):(c11cmn_combine(2)-c11cmn(4)+1);
                c11cmn_comTo12_lon=(c11cmn(1)-c11cmn_combine(1)+1):(c11cmn(3)-c11cmn_combine(1)+1);
            end

            lon_matl=c11cmn(1)-0.5:c11cmn(3)-0.5;
            lat_matl=c11cmn(2)+0.5:-1:c11cmn(4)+0.5;
            [lonlon_matl,latlat_matl]=meshgrid(lon_matl,lat_matl);
            %%
            % figure
            % font_Size=16;
            % % imagesc(RCFreDomi)
            % SHM=SHeatmap(RCFreDomi','Format','sq');
            % SHM=SHM.draw();
            % SHM.setText('FontSize',font_Size+2);
            %
            % for i=1:size(RCFreDomi',1)
            %     for j=1:size(RCFreDomi',2)
            %         if RC_noise(j,i)==1
            %             SHM.setTextMN(i,j,'String','x')
            %         else
            %             SHM.setTextMN(i,j,'String','')
            %         end
            %     end
            % end
            %
            % % SHM.setText('FontSize',font_Size+2);
            % colormap(mymap('coolwarm'))
            % colorbar


            %% plot
            %
            col_gap=[0 205 205]/255;col_est=[153 204 255]/255;col_out=[128 0 128]/255;
            col_art=[238 173 14]/255;
            % map_color=mymap("rainbow");
            font_Size=16;
            MarkSize=2;
            F_Position=[100,50,700,500];
            abc='abcdefghijklmnopqrstuvwxyz';

            % box plot for the energy distribution
            colorr_blind=[213,94,0;204,121,167;0,158,115;230,159,0;213,94,0;240,228,66;0,114,178;86,180,233]/255; % Blind-friendly colors

            ybou1=[0,1.19];
            F_gap=[.02 .03];F_marg_h=[.12 .08];F_marg_w=[.10 .10];

            used_eva=A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_eva;

            f1a=figure
            left_color=[0 0 0];
            right_color=[255 0 0]/255;
            set(f1a,'defaultAxesColorOrder',[left_color;right_color])
            ha=tight_subplot(1,1,F_gap,F_marg_h,F_marg_w);

            axes(ha(1));

            yyaxis left
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
            Y=[used_eva(1:S,[3,4,6]), used_eva(1:S,1)-used_eva(1:S,2)];

            bb=bar(X,Y,'stacked');
            set(bb(1),'FaceColor',colorr_blind(1,:));
            set(bb(2),'FaceColor',colorr_blind(3,:));
            % set(bb(3),'FaceColor',colorr_blind(2,:));
            set(bb(3),'FaceColor',colorr_blind(4,:));
            set(bb(4),'FaceColor',[169,169,169]/255);

            ylim(ybou1);
            xlim([0,S+1])

            ylabel('sum of normalized eigenvalues in M-SSA ($\sum\widetilde{\lambda}$)','Interpreter','latex')
            xlabel('Spherical slepian coefficients')

            yyaxis right

            y1=plot(1:S,V(1:S),'r','LineWidth',2);

            hold on
            % q1=quiver(S_shannon,CGJI_fill_MSSAcoffs_eva(S_shannon,1)+0.05, ...
            %     S_shannon,CGJI_fill_MSSAcoffs_eva(S_shannon,1),...
            %     'r','filled','linewidth',2)
            if S>=S_shannon
                text(S_shannon,used_eva(S_shannon,1)+0.05,'\downarrow', ...
                    'color','r','fontsize',22,'HorizontalAlignment','center')
                hold on
                text(S_shannon-0.5,used_eva(S_shannon,1)+0.1,'Shannon','fontsize',11)
            end

            ylim(ybou1);
            xlim([0,S+1])

            grid on

            ylabel('concentration ratio of spherical Slepian functions (\lambda)')

            %         legend({'long-term','annual','semi-annual','short-term','noise','\lambda'}, ...
            %             'NumColumns',3,'Position',[0.4,0.84,0.4,0.03],'box','off')
            legend({'long-term','annual','short-term','noise','\lambda'}, ...
                'NumColumns',3,'Position',[0.4,0.84,0.4,0.03],'box','off')

            title(['TP-S(' num2str(TP_ss) ') and TP-N(' num2str(TP_nn) ')'])

            set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.004,0.035]);

            xticks(get(gca, 'xtick'))
            set(gca, 'xticklabel', get(gca, 'xtick'))

            set(gcf,'Position',F_Position)

            tif_name1a=['MSSA_Eigenvalues_box_Frequency_V_TP' num2str(Turning_number) '_S' num2str(TP_ss) '_N' num2str(TP_nn) '_M' num2str(M_rec) '_Sig' num2str(Noise_SigLev) '.tif'];

            pause(1)
            print(f1a,'-dtiff','-r300',fullfile(fig_path_ALL,tif_name1a));

            %% Plot (also the noise frequency)
            %
            col_gap=[0 205 205]/255;col_est=[153 204 255]/255;col_out=[128 0 128]/255;
            col_art=[238 173 14]/255;
            % map_color=mymap("rainbow");
            font_Size=16;
            MarkSize=2;
            F_Position=[100,50,700,500];
            abc='abcdefghijklmnopqrstuvwxyz';

            % box plot for the energy distribution
            colorr_blind=[213,94,0;204,121,167;0,158,115;230,159,0;213,94,0;240,228,66;0,114,178;86,180,233]/255; % Blind-friendly colors

            ybou1=[0,1.19];
            F_gap=[.02 .03];F_marg_h=[.12 .08];F_marg_w=[.10 .10];

            used_eva=A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_eva;

            f1b=figure
            left_color=[0 0 0];
            right_color=[255 0 0]/255;
            set(f1b,'defaultAxesColorOrder',[left_color;right_color])
            ha=tight_subplot(1,1,F_gap,F_marg_h,F_marg_w);

            axes(ha(1));

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

            bb=bar(X,Y,'stacked');
            set(bb(1),'FaceColor',colorr_blind(1,:));
            set(bb(2),'FaceColor',colorr_blind(3,:));
            % set(bb(3),'FaceColor',colorr_blind(2,:));
            set(bb(3),'FaceColor',colorr_blind(4,:));

            hatchfill2(bb(4),'single','HatchAngle',60,'hatchcolor',colorr_blind(1,:),'HatchLineWidth',1.5);
            hatchfill2(bb(5),'single','HatchAngle',60,'hatchcolor',colorr_blind(3,:),'HatchLineWidth',1.5);
            hatchfill2(bb(6),'single','HatchAngle',60,'hatchcolor',colorr_blind(4,:),'HatchLineWidth',1.5);
            for b = 4:6
                %     bb(b).FaceColor = [0.8,0.8,0.8];
                bb(b).FaceColor = 'none';
            end


            % set(bb(4),'FaceColor',colorr_blind(1,:));
            % set(bb(5),'FaceColor',colorr_blind(3,:));
            % % set(bb(3),'FaceColor',colorr_blind(2,:));
            % set(bb(6),'FaceColor',colorr_blind(4,:));

            set(bb(7),'FaceColor',[169,169,169]/255);

            ylim(ybou1);
            xlim([0,S+1])

            grid on

            ylabel('sum of normalized eigenvalues in M-SSA ($\sum\widetilde{\lambda}$)','Interpreter','latex')
            xlabel('Spherical slepian coefficients')

            % yyaxis right

            % y1=plot(1:S,V(1:S),'r','LineWidth',2);

            % hold on
            % q1=quiver(S_shannon,CGJI_fill_MSSAcoffs_eva(S_shannon,1)+0.05, ...
            %     S_shannon,CGJI_fill_MSSAcoffs_eva(S_shannon,1),...
            %     'r','filled','linewidth',2)

            %         if S>=S_shannon
            %             text(S_shannon,used_eva(S_shannon,1)+0.05,'\downarrow', ...
            %                 'color','r','fontsize',22,'HorizontalAlignment','center')
            %             hold on
            %             text(S_shannon-0.5,used_eva(S_shannon,1)+0.1,'Shannon','fontsize',11)
            %         end

            % ylim(ybou1);

            % grid on

            % ylabel('concentration ratio of spherical Slepian functions (\lambda)')
            %
            %         legend([bb(1),bb(2),bb(3),bb(7)],{'long-term','annual','short-term','noise'}, ...
            %             'NumColumns',3,'Position',[0.4,0.84,0.4,0.03],'box','off')
            legend([bb(1),bb(2),bb(3),bb(7)],{'long-term','annual','short-term','residual'}, ...
                'NumColumns',3,'Location','northeast','box','off')
            tt=text(1,1.1,['\gamma = ' num2str((1-Noise_SigLev)*100) '%'],'FontSize',font_Size-2,...
                'verticalalignment','middle','horizontalalignment','left',...
                'backgroundColor','white','color','r','edgeColor','black');

            title({['Truncation number S(' num2str(TP_ss) ') and cutoff number N(' num2str(TP_nn) ')']},...
                'FontSize',font_Size-3)

            set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.004,0.035]);

            xticks(get(gca, 'xtick'))
            set(gca, 'xticklabel', get(gca, 'xtick'))

            set(gcf,'Position',F_Position)

            tif_name1b=['MSSA_Eigenvalues_box_noise_Frequency_V_TP' num2str(Turning_number) '_S' num2str(TP_ss) '_N' num2str(TP_nn) '_M' num2str(M_rec) '_Sig' num2str(Noise_SigLev) '.tif'];

            pause(1)
            print(f1b,'-dtiff','-r300',fullfile(fig_path_ALL,tif_name1b));

            %% select specific SSF, coefficient and some RCs (also show the SSF basis)
            % select_SSF=[3,7,15,16,21];
            %         if TP_ss==MaxNumChanges_V || TP_nn==MaxNumChanges_MSSA
            if TP_ss==MaxNumChanges_V && TP_nn==MaxNumChanges_MSSA
                select_S=1:S;
                % select_SSF=[3];

                Nc=6; % how many RCs do you want to show
                colorr20=mymap('tab20');
                for i=1:20
                    colorr_20(i,:)=colorr20(13*i-12,:);
                end
                col_gap=colorr_20(15,:);col_est=[153 204 255]/255;col_out=[128 0 128]/255;
                col_art=[238 173 14]/255;
                % map_color=mymap("rainbow");
                font_Size=16;
                MarkSize=2;
                F_Position=[2100,-250,1500,1100];
                abc='abcdefghijklmnopqrstuvwxyz';

                % box plot for the energy distribution
                % colorr_blind=[213,94,0;204,121,167;0,158,115;230,159,0;240,228,66;0,114,178;86,180,233]/255; % Blind-friendly colors
                colorr_blind=[0,114,178;204,121,167;0,158,115;230,159,0;213,94,0;240,228,66;86,180,233]/255; % Blind-friendly colors

                tt_use=floor(A5_Centers_Sort(1).use_months/12)+data_year_beg+mod(A5_Centers_Sort(1).use_months,12)/12-1/24;
                tt_fil=floor(fill_months/12)+data_year_beg+mod(fill_months,12)/12-1/24;
                tt_lea=floor(A5_Centers_Sort(1).missing_months/12)+data_year_beg+mod(A5_Centers_Sort(1).missing_months,12)/12-1/24;

                for i=1:numel(select_S)

                    Nc_plot=min(Nc,size(A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_RC{i},2));

                    f2=figure

                    set(gcf,'Position',F_Position)

                    SSF=select_S(i);

                    % plot the basis
                    F_gap=[.02 .03];F_marg_h=[.29 .24];F_marg_w=[.03 .78];

                    h0=tight_subplot(1,1,F_gap,F_marg_h,F_marg_w);

                    SSF_flat=reshape(squeeze(r_record(SSF,:,:)),[],1);

                    SSF_flat=normalize(SSF_flat,'norm',Inf);

                    SSF_nor=reshape(SSF_flat,size(r_example,1),size(r_example,2));

                    axes(h0(1));
                    %     m_proj('mercator','long',[c11cmn(1)-1.5, c11cmn(3)+0.5],'lat',[c11cmn(4)-0.5, c11cmn(2)+1.5]);
                    if strcmp(parrent_path(5:end),'GL')
                        m_proj('stereographic','long',(c11cmn(1)+c11cmn(3))/2,'lat',(c11cmn(4)+c11cmn(2))/2,...
                            'rad',20,'rec','on');
                        tt_lon=(c11cmn(3)-c11cmn(1))/5+c11cmn(1);tt_lat=c11cmn(4)+(c11cmn(2)-c11cmn(4))/12;
                    elseif strcmp(parrent_path(5:end),'Yang')
                        tt_lon=(c11cmn(3)-c11cmn(1))/20+c11cmn(1);tt_lat=c11cmn(4)+(c11cmn(2)-c11cmn(4))/6;
                        m_proj('mercator','long',[c11cmn(1)-0.5, c11cmn(3)-0.5],'lat',[c11cmn(4)+0.5, c11cmn(2)+0.5]);
                    else
                        tt_lon=(c11cmn(3)-c11cmn(1))/12+c11cmn(1);tt_lat=c11cmn(2)-(c11cmn(2)-c11cmn(4))/12;
                        m_proj('mercator','long',[c11cmn(1)-0.5, c11cmn(3)-0.5],'lat',[c11cmn(4)+0.5, c11cmn(2)+0.5]);
                    end
                    m_pcolor(lonlon_matl,latlat_matl,SSF_nor(c11cmn_comTo12_lat,c11cmn_comTo12_lon));

                    shading flat;
                    colormap(mymap("coolwarm"));
                    h=colorbar('h');
                    %     set(get(h,'xlabel'),'string','AVHRR SST Nov 1999');

                    m_coast
                    m_grid('box','fancy','tickdir','in');
                    %     m_grid('box','fancy','tickdir','in','xlabeldir', 'middle');
                    %     m_grid('box','fancy','tickdir','in','xlabeldir', 'middle');
                    m_line(XY_ori(:,1),XY_ori(:,2),'LineStyle','-','color','k','linewidth',1);
                    hold on
                    m_line(XY_buf(:,1),XY_buf(:,2),'LineStyle','--','color','r','linewidth',1);
                    % xlabel(['CC(' num2str(SSF) ')'])

                    hold on

                    m_text(tt_lon,tt_lat, ...
                        ['\lambda = ' num2str(A5_Centers_Sort(1).V(i),'%.2f')],'FontName','Times New Roman' ...
                        ,'fontsize',13,'color','b','verticalalignment','middle','horizontalalignment','left',...
                        'backgroundColor','white')

                    caxis([-1 1]);

                    if strcmp(parrent_path(5:end),'GL')
                        m_text(tt_lon,90, ...
                            ['\bf (a) Normalized slepian basis ($g_{' num2str(SSF) '}\left(\hat{\mathbf{r}}\right)$)'],...
                            'FontWeight','bold','FontSize',font_Size-1,'horizontalalignment','center','Interpreter','latex')
                    else
                        title(['\bf (a) Normalized slepian basis ($g_{' num2str(SSF) '}\left(\hat{\mathbf{r}}\right)$)'],...
                            'FontWeight','bold','FontSize',font_Size-1,'Interpreter','latex')
                    end

                    set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.004,0.035]);

                    % for the first original slepian coefficients
                    F_gap=[.02 .03];F_marg_h=[.82 .04];F_marg_w=[.3 .3];
                    ha=tight_subplot(1,1,F_gap,F_marg_h,F_marg_w);

                    ybou1=[-1,1];ybou2=[-20,20];

                    axes(ha(1));

                    pn=[];pn_name={};line_ins=[];
                    for ii=1:intit_num
                        line_ins(ii,:)=squeeze(A5_Centers_Sort(ii).fill_MSSAcoffs_original(:,i));
                        %             p1(ii)=plot(tt_fil,MSSAf_RC(:,i,ii),'color',map_colori(ii,:),'linewidth',2);
                        p1(ii)=plot(tt_fil,line_ins(ii,:),...
                            'color',colorr_blind(ii,:),'linewidth',2);
                        hold on

                        pn=[pn p1(ii)];pn_name{ii}=A5_Centers_Sort(ii).name(1:end-4);
                    end

                    Maxx=max(line_ins,[],'all');
                    Minn=min(line_ins,[],'all');

                    ybou1(1)=min(ybou1(1),Minn);
                    ybou1(2)=max(ybou1(2),Maxx);

                    Cap_posix1=tt_fil(5);Cap_posiy1=ybou1(2)-(ybou1(2)-ybou1(1))/10;
                    Cap_posix2=tt_fil(5);Cap_posiy2=ybou2(2)-(ybou2(2)-ybou2(1))/10;

                    for k=1:numel(tt_lea)
                        leak_tt_neibour=[tt_lea(k)-1/24,tt_lea(k),tt_lea(k)+1/24];
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

                    text(Cap_posix1,Cap_posiy1,['(' abc(2) ')'],'FontSize',font_Size-1);

                    %     if j==1
                    ylabel(['EWH (cm)'],'FontWeight','bold','FontSize',font_Size)
                    %     end

                    xlim([tt_fil(1)-0.5, tt_fil(end)+0.5])
                    ylim(ybou1)

                    pn_name{intit_num+1}='GRACE gap';
                    xticks(floor(tt_fil(1)):2:floor(tt_fil(end)))
                    xticklabels([]);

                    title(['\bf Spherical slepian coefficient ($d_{' num2str(SSF) '}$)'],'Interpreter','latex')

                    set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.004,0.035]);

                    %
                    F_gap=[.02 .03];F_marg_h=[.82 .04];F_marg_w=[.74 .02];
                    ha=tight_subplot(1,1,F_gap,F_marg_h,F_marg_w);

                    Fs = 12;                    % Sampling frequency
                    T = 1/Fs;                     % Sampling period
                    L = length(tt_fil);                     % Length of signal
                    t = (0:L-1)*T;                % Time vector
                    freq = Fs*(0:(L/2))/L;

                    ybou1=[0,0.5];ybou2=[0,1];

                    axes(ha(1));

                    Yss=zeros(size(line_ins(1,:)));
                    for ii=1:intit_num

                        Yss=Yss+squeeze(line_ins(ii,:));

                        Y = fft(squeeze(line_ins(ii,:)));
                        P2 = abs(Y/L);
                        P1 = P2(1:L/2+1);
                        P1(2:end-1) = 2*P1(2:end-1);

                        %             semilogx(1./freq,P1,'color',map_colori(ii,:),'LineWidth',1.2,'Marker','o','MarkerSize',MarkSize);
                        semilogx(1./freq,P1,'color',colorr_blind(ii,:),'LineWidth',1.2,'Marker','o','MarkerSize',MarkSize);
                        hold on

                        ybou1(2)=max(ybou1(2),max(P1));
                    end

                    xticks([3/12,6/12,1,4,10]);
                    xticklabels({'3','6','12','48','120'})
                    xticklabels([]);

                    ylim(ybou1)

                    Cap_posiy1=ybou1(2)-(ybou1(2)-ybou1(1))/10;
                    Cap_posiy2=ybou1(2)-(ybou1(2)-ybou1(1))/3;

                    text(0.12,Cap_posiy1,['(' abc(2+Nc_plot+1) ')'],'FontSize',font_Size-1);

                    title(['Power Spectrum'],'FontWeight','bold','FontSize',font_Size-1)

                    set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.007,0.035]);

                    F_gap=[.02 .03];F_marg_h=[.06 .22];F_marg_w=[.3 .3];
                    ha=tight_subplot(Nc_plot,1,F_gap,F_marg_h,F_marg_w);
                    map_colori=jet(4);

                    % Total_CGJI_EST_reconst_MSSA1=zeros(numel(tt_EST),intit_num);
                    for j=1:Nc_plot

                        ybou1=[-1,1];ybou2=[-20,20];

                        axes(ha(j));

                        pn=[];pn_name={};line_ins=[];
                        for ii=1:intit_num
                            line_ins(ii,:)=squeeze(A5_Centers_Sort(ii).fill_MSSAcoffs_RC{i}(:,j));
                            %             p1(ii)=plot(tt_fil,MSSAf_RC(:,i,ii),'color',map_colori(ii,:),'linewidth',2);
                            p1(ii)=plot(tt_fil,squeeze(A5_Centers_Sort(ii).fill_MSSAcoffs_RC{i}(:,j)),...
                                'color',colorr_blind(ii,:),'linewidth',2);
                            hold on

                            pn=[pn p1(ii)];pn_name{ii}=A5_Centers_Sort(ii).name(1:end-4);
                        end

                        Maxx=max(line_ins,[],'all');
                        Minn=min(line_ins,[],'all');

                        ybou1(1)=min(ybou1(1),Minn);
                        ybou1(2)=max(ybou1(2),Maxx);

                        Cap_posix1=tt_fil(5);Cap_posiy1=ybou1(2)-(ybou1(2)-ybou1(1))/7;
                        Cap_posix2=tt_fil(5);Cap_posiy2=ybou2(2)-(ybou2(2)-ybou2(1))/7;

                        for k=1:numel(tt_lea)
                            leak_tt_neibour=[tt_lea(k)-1/24,tt_lea(k),tt_lea(k)+1/24];
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

                        %                 text(Cap_posix1,Cap_posiy1,['(' abc(2+j) ')'],'FontSize',font_Size-1);
                        text(Cap_posix1,Cap_posiy1,['\bf (' abc(2+j) ') $\widetilde{\lambda}=' num2str(MSSAf_evalues_mm(j,i),'%.3f') '$'],...
                            'FontName','Times New Roman','FontSize',font_Size-1,'Interpreter','latex');

                        %     if j==1
                        ylabel(['RC-' num2str(j)],'FontWeight','bold','FontSize',font_Size)
                        %     end

                        xlim([tt_fil(1)-0.5, tt_fil(end)+0.5])
                        ylim(ybou1)

                        pn_name{intit_num+1}='GRACE gap';
                        xticks(floor(tt_fil(1)):2:floor(tt_fil(end)))

                        if j==1
                            if i==S_ol && S_ol>S
                                title(['\bf Leading reconstructed components (RCs) for inverted barometer'],'FontWeight','bold')
                            else
                                title(['\bf Leading reconstructed components (RCs) for $d_{' num2str(SSF) '}$'],'Interpreter','latex')
                            end
                            %             title(['M-SSA Reconstruction'],'FontWeight','bold','FontSize',font_Size-1)
                        end

                        if j==Nc_plot
                            xlabel('Time (year)')
                        else
                            xticklabels([]);
                        end

                        set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.004,0.035]);

                        %         title(['GRACE SSF MSSA-' num2str(i)],'FontWeight','bold')

                    end

                    F_gap=[.02 .03];F_marg_h=[.06 .22];F_marg_w=[.74 .02];
                    ha=tight_subplot(Nc_plot,1,F_gap,F_marg_h,F_marg_w);
                    map_colori=jet(4);

                    Fs = 12;                    % Sampling frequency
                    T = 1/Fs;                     % Sampling period
                    L = length(tt_fil);                     % Length of signal
                    t = (0:L-1)*T;                % Time vector
                    freq = Fs*(0:(L/2))/L;

                    RCTest=A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_RCTest{i};

                    ybou1=[0,0.5];ybou2=[0,1];
                    % Total_CGJI_EST_reconst_MSSA1=zeros(numel(tt_EST),intit_num);
                    for j=1:Nc_plot

                        axes(ha(j));

                        Yss=zeros(size(line_ins(1,:)));
                        for ii=1:intit_num
                            line_ins(ii,:)=squeeze(A5_Centers_Sort(ii).fill_MSSAcoffs_RC{i}(:,j));

                            Yss=Yss+squeeze(line_ins(ii,:));

                            Y = fft(squeeze(line_ins(ii,:)));
                            P2 = abs(Y/L);
                            P1 = P2(1:L/2+1);
                            P1(2:end-1) = 2*P1(2:end-1);

                            %             semilogx(1./freq,P1,'color',map_colori(ii,:),'LineWidth',1.2,'Marker','o','MarkerSize',MarkSize);
                            semilogx(1./freq,P1,'color',colorr_blind(ii,:),'LineWidth',1.2,'Marker','o','MarkerSize',MarkSize);
                            hold on

                            ybou1(2)=max(ybou1(2),max(P1));
                        end
                        %
                        %         hks(j) = kstest((Yss-mean(Yss))/std(Yss),'Alpha',0.05);
                        %         hli(j) = lillietest((Yss-mean(Yss))/std(Yss),'Alpha',0.05);

                        xticks([3/12,6/12,1,4,10]);
                        xticklabels({'3','6','12','48','120'})
                        if j<Nc_plot
                            xticklabels([]);
                        else
                            xlabel('Period (month)')
                        end

                        ybou1=get(gca,'ylim');
                        ylim(ybou1)

                        Cap_posiy1=ybou1(2)-(ybou1(2)-ybou1(1))/10;
                        Cap_posiy2=ybou1(2)-(ybou1(2)-ybou1(1))/3;

                        text(0.12,Cap_posiy1,['(' abc(2+j+Nc_plot+1) ')'],'FontSize',font_Size-1);

                        hcdf_fre=RCTest(4,:);
                        hks=RCTest(1,:);
                        hli=RCTest(2,:);
                        hcdf=RCTest(3,:);
                        pks=RCTest(8,:);
                        pli=RCTest(9,:);

                        test_x=2;
                        if hcdf_fre(j)<1
                            test_x=0.3;
                        end


                        %                 if hks(j)
                        %                     text(test_x,Cap_posiy1,['K-S: pass (5%)'],'FontSize',font_Size-3);
                        %                 else
                        %                     text(test_x,Cap_posiy1,['K-S: fail (5%)'],'FontSize',font_Size-3);
                        %                 end
                        %                 if hli(j)
                        %                     text(test_x,Cap_posiy2,['Lilliefors: pass (5%)'],'FontSize',font_Size-3);
                        %                 else
                        %                     text(test_x,Cap_posiy2,['Lilliefors: fail (5%)'],'FontSize',font_Size-3);
                        %                 end

                        if pks(j)<0.05
                            text(test_x,Cap_posiy1,['p(K-S)<0.05'],'FontSize',font_Size-3);
                        elseif pks(j)<0.1
                            text(test_x,Cap_posiy1,['p(K-S)<0.1'],'FontSize',font_Size-3);
                        elseif pks(j)<0.3
                            text(test_x,Cap_posiy1,['p(K-S)<0.3'],'FontSize',font_Size-3,'color',[0.4940 0.1840 0.5560]);
                        else
                            text(test_x,Cap_posiy1,['p(K-S)>0.3'],'FontSize',font_Size-3,'color','r');
                        end

                        if pli(j)<0.05
                            text(test_x,Cap_posiy2,['p(Lilliefors)<0.05'],'FontSize',font_Size-3);
                        elseif pli(j)<0.1
                            text(test_x,Cap_posiy2,['p(Lilliefors)<0.1'],'FontSize',font_Size-3);
                        elseif pli(j)<0.3
                            text(test_x,Cap_posiy2,['p(Lilliefors)<0.3'],'FontSize',font_Size-3,'color',[0.4940 0.1840 0.5560]);
                        else
                            text(test_x,Cap_posiy2,['p(Lilliefors)>0.3'],'FontSize',font_Size-3,'color','r');
                        end

                        set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.007,0.035]);

                    end

                    if any(A5_Centers_Sort(1).missing_months)
                        pn_name{intit_num+1}='GRACE gap';
                        l1=legend([pn,he_gap],pn_name,'NumColumns',2,'box','off','Location','northwest')
                    else
                        l1=legend([pn],pn_name,'NumColumns',2,'box','off','Location','northwest')
                    end
                    set(l1,'Position',[0.01,0.2,0.2,0.01],'box','off','fontsize',font_Size+1)

                    tif_name2=['MSSA_Decompose_CC' num2str(SSF) '.tif'];
                    if i==4
                        print(f2,'-dtiff','-r500',fullfile(fig_path_ALL,tif_name2));
                    else
                        print(f2,'-dtiff','-r300',fullfile(fig_path_ALL,tif_name2));
                    end

                end

            end
            %%


            %% Reconstruct the gridded EWH (MSL)
            % Some criteria for the calculation is that:
            % we have the overall Mass and the gridded EWH(OBP) through the slepian
            % calculation, then we need to calculate:
            % (1) basin-averaged EWH (or basin-averaged MSL): from overall Mass
            % (2) gridded Mass (or gridded MSL): from gridded EWH

            % In addition, there is another process (M-SSA). Then, we will first
            % calculate the corresponding overall MASS and the gridded EWH(OBP), then
            % other variables as mentioned above.

            % Besides, for the original SHC product, only the overall Mass and the
            % gridded EWH(OBP) will be provided.

            fitwhat=[3 365.25 365.25/2]; %linear trend and an annual and a semiannual harmonic term

            [~,N_SHC]=size(A5_Centers_Sort(1).fill_SHCfilldelta);
            [in, on]=check_polygon_in(XY,lonlon,latlat);

            % area-weighted should be used to more precisely calculate the regional
            % basin-averaged variables
            for i=1:length(lat)
                for j=1:length(lon)
                    c11cmn_area(i,j)=areaquad(lat(i)-0.5,lon(j)-0.5,lat(i)+0.5,lon(j)+0.5);
                end
            end

            %% M-SSA (Slepian) part
            % we first calculate the overall Mass and the gridded EWH(OBP) through the
            % slepian calculation, then we need to calculate:
            % (1) basin-averaged EWH (or basin-averaged MSL): from overall Mass
            % (2) gridded Mass (or gridded MSL): from gridded EWH

            % do selection criteria for both K-S and lillietest tests to rejected noise
            % this is for the basin-averaged Mass or EWH (MSL)
            for ins=1:intit_num+1

                % the original coefficients without MSSA has already been calculated
                % for the spatial and temporal results.
                fill_deltacoffs_ins=A5_Centers_Sort(ins).fill_deltacoffs;
                fill_signaldelta_ins=A5_Centers_Sort(ins).fill_signaldelta;

                % (spatial)
                % the gap-filling part with MSSA
                fill_MSSAcoffs_ins=A5_Centers_Sort(ins).fill_MSSAcoffs;
                fill_MSSAcoffs_both_ins=A5_Centers_Sort(ins).fill_MSSAcoffs_both;

                for i=1:numel(fill_months)
                    sp_ewh_mssa_reconst=0;sp_ewh_mssa_both_reconst=0;
                    sp_ewh_mssa_resid_reconst=0;sp_ewh_mssa_bothfail_reconst=0;
                    for j=1:S
                        %             [r,lon,lat,Plm,degres]=plm2xyz(CC{j},1,c11cmn,60);
                        r=squeeze(r_record(j,:,:));

                        sp_ewh_mssa_reconst=sp_ewh_mssa_reconst+r*fill_MSSAcoffs_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm for sea water (/ 1000 kg/m*3 * 100)
                        sp_ewh_mssa_both_reconst=sp_ewh_mssa_both_reconst+r*fill_MSSAcoffs_both_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm for sea water (/ 1000 kg/m*3 * 100)
                        sp_ewh_mssa_resid_reconst=sp_ewh_mssa_resid_reconst+r*(fill_deltacoffs_ins(i,j)'-fill_MSSAcoffs_ins(i,j)') ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm for sea water (/ 1000 kg/m*3 * 100)
                        sp_ewh_mssa_bothfail_reconst=sp_ewh_mssa_bothfail_reconst+r*(fill_MSSAcoffs_ins(i,j)'-fill_MSSAcoffs_both_ins(i,j)') ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm for sea water (/ 1000 kg/m*3 * 100)Grid_EWH

                        % notice that regional EWH is the area-weighted mean, while
                        % regional Mass is the sum of all mass
                        A5_Centers_Sort(ins).fill_uptoS_EWH_MSSA(j,i)=sum(sp_ewh_mssa_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                        A5_Centers_Sort(ins).fill_uptoS_EWH_MSSA_both(j,i)=sum(sp_ewh_mssa_both_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                        A5_Centers_Sort(ins).fill_uptoS_EWH_MSSA_resid(j,i)=sum(sp_ewh_mssa_resid_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                        A5_Centers_Sort(ins).fill_uptoS_EWH_MSSA_bothfail(j,i)=sum(sp_ewh_mssa_bothfail_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));

                    end

                    A5_Centers_Sort(ins).fill_Grid_EWH_MSSA(i,:,:)=sp_ewh_mssa_reconst;
                    A5_Centers_Sort(ins).fill_Grid_EWH_MSSA_both(i,:,:)=sp_ewh_mssa_both_reconst;
                    A5_Centers_Sort(ins).fill_Grid_EWH_MSSA_resid(i,:,:)=sp_ewh_mssa_resid_reconst;
                    A5_Centers_Sort(ins).fill_Grid_EWH_MSSA_bothfail(i,:,:)=sp_ewh_mssa_bothfail_reconst;
                    A5_Centers_Sort(ins).fill_Grid_MASS_MSSA(i,:,:)=sp_ewh_mssa_reconst/100.*c11cmn_area*4*pi*6370000^2*10^3/10^3/10^9; % from cm to Gt
                    A5_Centers_Sort(ins).fill_Grid_MASS_MSSA_both(i,:,:)=sp_ewh_mssa_both_reconst/100.*c11cmn_area*4*pi*6370000^2*10^3/10^3/10^9; % from cm to Gt
                    A5_Centers_Sort(ins).fill_Grid_MASS_MSSA_resid(i,:,:)=sp_ewh_mssa_resid_reconst/100.*c11cmn_area*4*pi*6370000^2*10^3/10^3/10^9; % from cm to Gt
                    A5_Centers_Sort(ins).fill_Grid_MASS_MSSA_bothfail(i,:,:)=sp_ewh_mssa_bothfail_reconst/100.*c11cmn_area*4*pi*6370000^2*10^3/10^3/10^9; % from cm to Gt

                    % notice that regional EWH is the area-weighted mean, while
                    % regional Mass is the sum of all mass
                    A5_Centers_Sort(ins).fill_EWH_MSSA(i)=sum(sp_ewh_mssa_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                    A5_Centers_Sort(ins).fill_EWH_MSSA_both(i)=sum(sp_ewh_mssa_both_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                    A5_Centers_Sort(ins).fill_EWH_MSSA_resid(i)=sum(sp_ewh_mssa_resid_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                    A5_Centers_Sort(ins).fill_EWH_MSSA_bothfail(i)=sum(sp_ewh_mssa_bothfail_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                    A5_Centers_Sort(ins).fill_MASS_MSSA(i)=A5_Centers_Sort(ins).fill_EWH_MSSA(i)/100*BasinArea*10^3/10^3/10^9; % change from cm to Gt
                    A5_Centers_Sort(ins).fill_MASS_MSSA_both(i)=A5_Centers_Sort(ins).fill_EWH_MSSA_both(i)/100*BasinArea*10^3/10^3/10^9; % change from cm to Gt
                    A5_Centers_Sort(ins).fill_MASS_MSSA_resid(i)=A5_Centers_Sort(ins).fill_EWH_MSSA_resid(i)/100*BasinArea*10^3/10^3/10^9; % change from cm to Gt
                    A5_Centers_Sort(ins).fill_MASS_MSSA_bothfail(i)=A5_Centers_Sort(ins).fill_EWH_MSSA_bothfail(i)/100*BasinArea*10^3/10^3/10^9; % change from cm to Gt
                end

                % this is to calculate the signal (fit component) and residual (unfit
                % component) to be comsistent with previous process on slepian
                [fill_MSSA_signaldelta_ins,fill_MSSA_residdelta_ins,~,~]=slept2resid_m(fill_MSSAcoffs_ins,...
                    fill_dates,fitwhat,[],[],[],[],[],[]); % mm (values equivalent to kg/m^2)
                [fill_MSSA_bothsignaldelta_ins,fill_MSSA_bothresiddelta_ins,~,~]=slept2resid_m(fill_MSSAcoffs_both_ins,...
                    fill_dates,fitwhat,[],[],[],[],[],[]); % mm (values equivalent to kg/m^2)

                for i=1:numel(fill_months)
                    sp_ewhsig_mssa_reconst=0;sp_ewhres_mssa_reconst=0;
                    sp_ewhsig_mssa_both_reconst=0;sp_ewhres_mssa_both_reconst=0;
                    for j=1:S
                        %             [r,lon,lat,Plm,degres]=plm2xyz(CC{j},1,c11cmn,60);
                        r=squeeze(r_record(j,:,:));

                        sp_ewhsig_mssa_reconst=sp_ewhsig_mssa_reconst+r*fill_MSSA_signaldelta_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm for sea water (/ 1000 kg/m*3 * 100)
                        sp_ewhres_mssa_reconst=sp_ewhres_mssa_reconst+r*fill_MSSA_residdelta_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm for sea water (/ 1000 kg/m*3 * 100)
                        sp_ewhsig_mssa_both_reconst=sp_ewhsig_mssa_both_reconst+r*fill_MSSA_bothsignaldelta_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm for sea water (/ 1000 kg/m*3 * 100)
                        sp_ewhres_mssa_both_reconst=sp_ewhres_mssa_both_reconst+r*fill_MSSA_bothresiddelta_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm for sea water (/ 1000 kg/m*3 * 100)

                        % notice that regional EWH is the area-weighted mean, while
                        % regional Mass is the sum of all mass
                        A5_Centers_Sort(ins).fill_uptoS_EWHsig_MSSA(j,i)=sum(sp_ewhsig_mssa_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                        A5_Centers_Sort(ins).fill_uptoS_EWHres_MSSA(j,i)=sum(sp_ewhres_mssa_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                        A5_Centers_Sort(ins).fill_uptoS_EWHsig_MSSA_both(j,i)=sum(sp_ewhsig_mssa_both_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                        A5_Centers_Sort(ins).fill_uptoS_EWHres_MSSA_both(j,i)=sum(sp_ewhres_mssa_both_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));

                    end

                    A5_Centers_Sort(ins).fill_Grid_EWHsig_MSSA(i,:,:)=sp_ewhsig_mssa_reconst;
                    A5_Centers_Sort(ins).fill_Grid_EWHres_MSSA(i,:,:)=sp_ewhres_mssa_reconst;
                    A5_Centers_Sort(ins).fill_Grid_EWHsig_MSSA_both(i,:,:)=sp_ewhsig_mssa_both_reconst;
                    A5_Centers_Sort(ins).fill_Grid_EWHres_MSSA_both(i,:,:)=sp_ewhres_mssa_both_reconst;
                    A5_Centers_Sort(ins).fill_Grid_MASSsig_MSSA(i,:,:)=sp_ewhsig_mssa_reconst/100.*c11cmn_area*4*pi*6370000^2*10^3/10^3/10^9; % from cm to Gt
                    A5_Centers_Sort(ins).fill_Grid_MASSres_MSSA(i,:,:)=sp_ewhres_mssa_reconst/100.*c11cmn_area*4*pi*6370000^2*10^3/10^3/10^9; % from cm to Gt
                    A5_Centers_Sort(ins).fill_Grid_MASSsig_MSSA_both(i,:,:)=sp_ewhsig_mssa_both_reconst/100.*c11cmn_area*4*pi*6370000^2*10^3/10^3/10^9; % from cm to Gt
                    A5_Centers_Sort(ins).fill_Grid_MASSres_MSSA_both(i,:,:)=sp_ewhres_mssa_both_reconst/100.*c11cmn_area*4*pi*6370000^2*10^3/10^3/10^9; % from cm to Gt

                    % notice that regional EWH is the area-weighted mean, while
                    % regional Mass is the sum of all mass
                    A5_Centers_Sort(ins).fill_EWHsig_MSSA(i)=sum(sp_ewhsig_mssa_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                    A5_Centers_Sort(ins).fill_EWHres_MSSA(i)=sum(sp_ewhres_mssa_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                    A5_Centers_Sort(ins).fill_EWHsig_MSSA_both(i)=sum(sp_ewhsig_mssa_both_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                    A5_Centers_Sort(ins).fill_EWHres_MSSA_both(i)=sum(sp_ewhres_mssa_both_reconst(in).*c11cmn_area(in))/sum(c11cmn_area(in));
                    A5_Centers_Sort(ins).fill_MASSsig_MSSA(i)=A5_Centers_Sort(ins).fill_EWHsig_MSSA(i)/100*BasinArea*10^3/10^3/10^9; % change from cm to Gt
                    A5_Centers_Sort(ins).fill_MASSres_MSSA(i)=A5_Centers_Sort(ins).fill_EWHres_MSSA(i)/100*BasinArea*10^3/10^3/10^9; % change from cm to Gt
                    A5_Centers_Sort(ins).fill_MASSsig_MSSA_both(i)=A5_Centers_Sort(ins).fill_EWHsig_MSSA_both(i)/100*BasinArea*10^3/10^3/10^9; % change from cm to Gt
                    A5_Centers_Sort(ins).fill_MASSres_MSSA_both(i)=A5_Centers_Sort(ins).fill_EWHres_MSSA_both(i)/100*BasinArea*10^3/10^3/10^9; % change from cm to Gt
                end

                if S_ol>S

                    for i=1:numel(fill_months)

                        A5_Centers_Sort(ins).fill_MSL_MSSA(i)=A5_Centers_Sort(ins).fill_EWH_MSSA(i)*1000/1028-...   % sea water
                            fill_MSSAcoffs_ins(i,S_ol)*1000/1028/10; % IB (unit:mmH2O, kg/m^2) to cm
                        A5_Centers_Sort(ins).fill_MSL_MSSA_both(i)=A5_Centers_Sort(ins).fill_EWH_MSSA_both(i)*1000/1028-...   % sea water
                            fill_MSSAcoffs_both_ins(i,S_ol)*1000/1028/10; % IB (unit:mmH2O, kg/m^2) to cm
                        A5_Centers_Sort(ins).fill_MSL_MSSA_resid(i)=A5_Centers_Sort(ins).fill_EWH_MSSA_resid(i)*1000/1028-...   % sea water
                            (A5_Centers_Sort(ins).fill_IBdeltacoffs(i,1)-fill_MSSAcoffs_ins(i,S_ol))*1000/1028/10; % IB (unit:mmH2O, kg/m^2) to cm
                        A5_Centers_Sort(ins).fill_MSL_MSSA_bothfail(i)=A5_Centers_Sort(ins).fill_EWH_MSSA_bothfail(i)*1000/1028-...   % sea water
                            (fill_MSSAcoffs_ins(i,S_ol)-fill_MSSAcoffs_both_ins(i,S_ol))*1000/1028/10; % IB (unit:mmH2O, kg/m^2) to cm

                        A5_Centers_Sort(ins).fill_MSLsig_MSSA(i)=A5_Centers_Sort(ins).fill_MASSsig_MSSA(i)*10^3*10^9/10^3/BasinArea*100-... %from Gt to cm
                            fill_MSSA_signaldelta_ins(i,S_ol)*1000/1028/10;%from Gt to cm
                        A5_Centers_Sort(ins).fill_MSLres_MSSA(i)=A5_Centers_Sort(ins).fill_MASSres_MSSA(i)*10^3*10^9/10^3/BasinArea*100-... %from Gt to cm
                            fill_MSSA_residdelta_ins(i,S_ol)*1000/1028/10;%from Gt to cm
                        A5_Centers_Sort(ins).fill_MSLsig_MSSA_both(i)=A5_Centers_Sort(ins).fill_MASSsig_MSSA_both(i)*10^3*10^9/10^3/BasinArea*100-... %from Gt to cm
                            fill_MSSA_bothsignaldelta_ins(i,S_ol)*1000/1028/10;%from Gt to cm
                        A5_Centers_Sort(ins).fill_MSLres_MSSA_both(i)=A5_Centers_Sort(ins).fill_MASSres_MSSA_both(i)*10^3*10^9/10^3/BasinArea*100-... %from Gt to cm
                            fill_MSSA_bothresiddelta_ins(i,S_ol)*1000/1028/10;%from Gt to cm

                        % transfer EWH to OBP and then substract IB to gain MSL
                        A5_Centers_Sort(ins).fill_Grid_MSL_MSSA(i,:,:)=A5_Centers_Sort(ins).fill_Grid_EWH_MSSA(i,:,:)*1000/1028-... % sea water
                            fill_MSSAcoffs_ins(i,S_ol)*1000/1028/10; % IB (unit:mmH2O, kg/m^2) to cm
                        A5_Centers_Sort(ins).fill_Grid_MSL_MSSA_both(i,:,:)=A5_Centers_Sort(ins).fill_Grid_EWH_MSSA_both(i,:,:)*1000/1028-... % sea water
                            fill_MSSAcoffs_both_ins(i,S_ol)*1000/1028/10; % IB (unit:mmH2O, kg/m^2)
                        A5_Centers_Sort(ins).fill_Grid_MSL_MSSA_resid(i,:,:)=A5_Centers_Sort(ins).fill_Grid_EWH_MSSA_resid(i,:,:)*1000/1028-... % sea water
                            (A5_Centers_Sort(ins).fill_IBdeltacoffs(i,1)-fill_MSSAcoffs_ins(i,S_ol))*1000/1028/10; % IB (unit:mmH2O, kg/m^2)
                        A5_Centers_Sort(ins).fill_Grid_MSL_MSSA_bothfail(i,:,:)=A5_Centers_Sort(ins).fill_Grid_EWH_MSSA_bothfail(i,:,:)*1000/1028-... % sea water
                            (fill_MSSAcoffs_ins(i,S_ol)-fill_MSSAcoffs_both_ins(i,S_ol))*1000/1028/10; % IB (unit:mmH2O, kg/m^2)

                        A5_Centers_Sort(ins).fill_Grid_MSLsig_MSSA_both(i,:,:)=A5_Centers_Sort(ins).fill_Grid_EWHsig_MSSA_both(i,:,:)*1000/1028-... % sea water
                            fill_MSSA_bothsignaldelta_ins(i,S_ol)*1000/1028/10; % IB (unit:mmH2O, kg/m^2) to cm
                        A5_Centers_Sort(ins).fill_Grid_MSLres_MSSA_both(i,:,:)=A5_Centers_Sort(ins).fill_Grid_EWHres_MSSA_both(i,:,:)*1000/1028-... % sea water
                            fill_MSSA_bothresiddelta_ins(i,S_ol)*1000/1028/10; % IB (unit:mmH2O, kg/m^2)
                    end

                end

            end

            %% STD and RMS

            % calculate the STD for the reject and residuals (maybe?)
            % notice that here we do not consider IB, and only consider the EWH (e.g., use 1000 kg/m^3)
            A5_Centers_STD=[];A5_Centers_RMS=[];
            for ins=1:intit_num+1

                % the original coefficients without MSSA has already been calculated
                % for the spatial and temporal results.
                fill_deltacoffs_ins=A5_Centers_Sort(ins).fill_deltacoffs;
                fill_signaldelta_ins=A5_Centers_Sort(ins).fill_signaldelta;

                % (spatial)
                % the gap-filling part with MSSA
                fill_MSSAcoffs_ins=A5_Centers_Sort(ins).fill_MSSAcoffs;
                fill_MSSAcoffs_both_ins=A5_Centers_Sort(ins).fill_MSSAcoffs_both;

                % this is to calculate the signal (fit component) and residual (unfit
                % component) to be comsistent with previous process on slepian
                [fill_MSSA_signaldelta_ins,fill_MSSA_residdelta_ins,~,~]=slept2resid_m(fill_MSSAcoffs_ins,...
                    fill_dates,fitwhat,[],[],[],[],[],[]); % mm (values equivalent to kg/m^2)
                [fill_MSSA_bothsignaldelta_ins,fill_MSSA_bothresiddelta_ins,~,~]=slept2resid_m(fill_MSSAcoffs_both_ins,...
                    fill_dates,fitwhat,[],[],[],[],[],[]); % mm (values equivalent to kg/m^2)


                fill_Grid_EWH_STD=zeros(S,size(r_example,1),size(r_example,2));
                fill_Grid_EWHsig_STD=zeros(S,size(r_example,1),size(r_example,2));
                fill_Grid_EWHres_STD=zeros(S,size(r_example,1),size(r_example,2));
                fill_Grid_EWH_uptoS_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWHsig_uptoS_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWHres_uptoS_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));

                fill_Grid_EWH_MSSA_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWH_MSSA_both_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWHsig_MSSA_both_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWHres_MSSA_both_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWH_MSSA_bothfail_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWH_MSSA_resid_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));

                fill_Grid_EWH_MSSA_uptoS_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWH_MSSA_both_uptoS_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWHsig_MSSA_both_uptoS_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWHres_MSSA_both_uptoS_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWH_MSSA_bothfail_uptoS_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWH_MSSA_resid_uptoS_STD=zeros(numel(fill_months),size(r_example,1),size(r_example,2));

                fill_Grid_EWH_RMS=zeros(S,size(r_example,1),size(r_example,2));
                fill_Grid_EWHsig_RMS=zeros(S,size(r_example,1),size(r_example,2));
                fill_Grid_EWHres_RMS=zeros(S,size(r_example,1),size(r_example,2));
                fill_Grid_EWH_uptoS_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWHsig_uptoS_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWHres_uptoS_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));

                fill_Grid_EWH_MSSA_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWH_MSSA_both_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWHsig_MSSA_both_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWHres_MSSA_both_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWH_MSSA_bothfail_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWH_MSSA_resid_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));

                fill_Grid_EWH_MSSA_uptoS_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWH_MSSA_both_uptoS_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWHsig_MSSA_both_uptoS_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWHres_MSSA_both_uptoS_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWH_MSSA_bothfail_uptoS_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                fill_Grid_EWH_MSSA_resid_uptoS_RMS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));

                sp_ewh_uptoS=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                sp_ewh_uptoS_signal=zeros(numel(fill_months),size(r_example,1),size(r_example,2));

                sp_ewh_uptoS_mssa=zeros(numel(fill_months),size(r_example,1),size(r_example,2));

                sp_ewh_uptoS_mssa_both=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                sp_ewh_uptoS_mssa_both_signal=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                sp_ewh_uptoS_mssa_both_resid=zeros(numel(fill_months),size(r_example,1),size(r_example,2));

                for j=1:S
                    sp_ewh=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                    sp_ewh_signal=zeros(numel(fill_months),size(r_example,1),size(r_example,2));

                    sp_ewh_mssa=zeros(numel(fill_months),size(r_example,1),size(r_example,2));

                    sp_ewh_mssa_both=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                    sp_ewh_mssa_both_signal=zeros(numel(fill_months),size(r_example,1),size(r_example,2));
                    sp_ewh_mssa_both_resid=zeros(numel(fill_months),size(r_example,1),size(r_example,2));

                    r=squeeze(r_record(j,:,:));

                    for i=1:numel(fill_months)
                        sp_ewh(i,:,:)=r*fill_deltacoffs_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
                        sp_ewh_signal(i,:,:)=r*fill_signaldelta_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
                        sp_ewh_uptoS(i,:,:)=squeeze(sp_ewh_uptoS(i,:,:))+r*fill_deltacoffs_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
                        sp_ewh_uptoS_signal(i,:,:)=squeeze(sp_ewh_uptoS_signal(i,:,:))+r*fill_signaldelta_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)

                        sp_ewh_mssa(i,:,:)=r*fill_MSSAcoffs_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
                        sp_ewh_mssa_both(i,:,:)=r*fill_MSSAcoffs_both_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
                        sp_ewh_mssa_both_signal(i,:,:)=r*fill_MSSA_bothsignaldelta_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
                        sp_ewh_mssa_both_resid(i,:,:)=r*fill_MSSA_bothresiddelta_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)

                        sp_ewh_uptoS_mssa(i,:,:)=squeeze(sp_ewh_uptoS_mssa(i,:,:))+r*fill_MSSAcoffs_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
                        sp_ewh_uptoS_mssa_both(i,:,:)=squeeze(sp_ewh_uptoS_mssa_both(i,:,:))+r*fill_MSSAcoffs_both_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
                        sp_ewh_uptoS_mssa_both_signal(i,:,:)=squeeze(sp_ewh_uptoS_mssa_both_signal(i,:,:))+r*fill_MSSA_bothsignaldelta_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
                        sp_ewh_uptoS_mssa_both_resid(i,:,:)=squeeze(sp_ewh_uptoS_mssa_both_resid(i,:,:))+r*fill_MSSA_bothresiddelta_ins(i,j)' ... % (you can manually adjust this period)
                            /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
                    end

                    for p=1:size(lonlon,1)
                        for q=1:size(lonlon,2)
                            A5_Centers_STD(ins).fill_Grid_EWH_STD(j,p,q)=std(squeeze(sp_ewh(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWHsig_STD(j,p,q)=std(squeeze(sp_ewh_signal(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWHres_STD(j,p,q)=std(squeeze(sp_ewh(:,p,q))-squeeze(sp_ewh_signal(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWH_uptoS_STD(j,p,q)=std(squeeze(sp_ewh_uptoS(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWHsig_uptoS_STD(j,p,q)=std(squeeze(sp_ewh_uptoS_signal(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWHres_uptoS_STD(j,p,q)=std(squeeze(sp_ewh_uptoS(:,p,q))-squeeze(sp_ewh_uptoS_signal(:,p,q)));

                            A5_Centers_STD(ins).fill_Grid_EWH_MSSA_STD(j,p,q)=std(squeeze(sp_ewh_mssa(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWH_MSSA_uptoS_STD(j,p,q)=std(squeeze(sp_ewh_uptoS_mssa(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWH_MSSA_both_STD(j,p,q)=std(squeeze(sp_ewh_mssa_both(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWH_MSSA_both_uptoS_STD(j,p,q)=std(squeeze(sp_ewh_uptoS_mssa_both(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWHsig_MSSA_both_STD(j,p,q)=std(squeeze(sp_ewh_mssa_both_signal(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWHsig_MSSA_both_uptoS_STD(j,p,q)=std(squeeze(sp_ewh_uptoS_mssa_both_signal(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWHres_MSSA_both_STD(j,p,q)=std(squeeze(sp_ewh_mssa_both_resid(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWHres_MSSA_both_uptoS_STD(j,p,q)=std(squeeze(sp_ewh_uptoS_mssa_both_resid(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWH_MSSA_bothfail_STD(j,p,q)=std(squeeze(sp_ewh_mssa(:,p,q)-sp_ewh_mssa_both(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWH_MSSA_bothfail_uptoS_STD(j,p,q)=std(squeeze(sp_ewh_uptoS_mssa(:,p,q)-sp_ewh_uptoS_mssa_both(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWH_MSSA_resid_STD(j,p,q)=std(squeeze(sp_ewh(:,p,q)-sp_ewh_mssa(:,p,q)));
                            A5_Centers_STD(ins).fill_Grid_EWH_MSSA_resid_uptoS_STD(j,p,q)=std(squeeze(sp_ewh_uptoS(:,p,q)-sp_ewh_uptoS_mssa(:,p,q)));

                            A5_Centers_RMS(ins).fill_Grid_EWH_RMS(j,p,q)=rms(squeeze(sp_ewh(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWHsig_RMS(j,p,q)=rms(squeeze(sp_ewh_signal(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWHres_RMS(j,p,q)=rms(squeeze(sp_ewh(:,p,q))-squeeze(sp_ewh_signal(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWH_uptoS_RMS(j,p,q)=rms(squeeze(sp_ewh_uptoS(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWHsig_uptoS_RMS(j,p,q)=rms(squeeze(sp_ewh_uptoS_signal(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWHres_uptoS_RMS(j,p,q)=rms(squeeze(sp_ewh_uptoS(:,p,q))-squeeze(sp_ewh_uptoS_signal(:,p,q)));

                            A5_Centers_RMS(ins).fill_Grid_EWH_MSSA_RMS(j,p,q)=rms(squeeze(sp_ewh_mssa(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWH_MSSA_uptoS_RMS(j,p,q)=rms(squeeze(sp_ewh_uptoS_mssa(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWH_MSSA_both_RMS(j,p,q)=rms(squeeze(sp_ewh_mssa_both(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWH_MSSA_both_uptoS_RMS(j,p,q)=rms(squeeze(sp_ewh_uptoS_mssa_both(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWHsig_MSSA_both_RMS(j,p,q)=rms(squeeze(sp_ewh_mssa_both_signal(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWHsig_MSSA_both_uptoS_RMS(j,p,q)=rms(squeeze(sp_ewh_uptoS_mssa_both_signal(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWHres_MSSA_both_RMS(j,p,q)=rms(squeeze(sp_ewh_mssa_both_resid(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWHres_MSSA_both_uptoS_RMS(j,p,q)=rms(squeeze(sp_ewh_uptoS_mssa_both_resid(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWH_MSSA_bothfail_RMS(j,p,q)=rms(squeeze(sp_ewh_mssa(:,p,q)-sp_ewh_mssa_both(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWH_MSSA_bothfail_uptoS_RMS(j,p,q)=rms(squeeze(sp_ewh_uptoS_mssa(:,p,q)-sp_ewh_uptoS_mssa_both(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWH_MSSA_resid_RMS(j,p,q)=rms(squeeze(sp_ewh(:,p,q)-sp_ewh_mssa(:,p,q)));
                            A5_Centers_RMS(ins).fill_Grid_EWH_MSSA_resid_uptoS_RMS(j,p,q)=rms(squeeze(sp_ewh_uptoS(:,p,q)-sp_ewh_uptoS_mssa(:,p,q)));
                        end
                    end

                end

            end
            %

            %% Plot the temporal final mass change (Gt) or EWH (or MSL) change (cm)

            col_gap=[0 205 205]/255;col_est=[153 204 255]/255;col_out=[128 0 128]/255;
            col_art=[238 173 14]/255;
            % map_color=mymap("rainbow");
            font_Size=16;
            MarkSize=2;
            F_Position=[100,50,1500,900];
            abc='abcdefghijklmnopqrstuvwxyz';

            % box plot for the energy distribution
            % colorr_blind=[213,94,0;204,121,167;0,158,115;230,159,0;240,228,66;0,114,178;86,180,233]/255; % Blind-friendly colors
            colorr_blind=[0,114,178;204,121,167;0,158,115;230,159,0;213,94,0;240,228,66;86,180,233]/255; % Blind-friendly colors

            tt_use=floor(A5_Centers_Sort(1).use_months/12)+data_year_beg+mod(A5_Centers_Sort(1).use_months,12)/12-1/24;
            tt_fil=floor(fill_months/12)+data_year_beg+mod(fill_months,12)/12-1/24;
            tt_lea=floor(A5_Centers_Sort(1).missing_months/12)+data_year_beg+mod(A5_Centers_Sort(1).missing_months,12)/12-1/24;

            for ins=1:intit_num+1
                % the SHC mass
                %             [fit_MASS,delta_MASS,totalparams_MASS,paramerrors_MASS] = timeseriesfit([fill_dates' A5_Centers_Sort(ins).fill_MASS_SHC'], ...
                %                 var(A5_Centers_Sort(ins).fill_MASSres_SHC),1,1);
                %             A5_Centers_Fit(ins).fill_MASS_SHC_slope=totalparams_MASS(2)*365.25;
                %             A5_Centers_Fit(ins).fill_MASS_SHC_fit = [fill_dates' fit_MASS delta_MASS];
                %             A5_Centers_Fit(ins).fill_MASS_SHC_paramerrors = paramerrors_MASS(2)*365.25; % this (365.25) is special for 'timeseriesfit' results.

                % the Slepian mass
                %             [fit_MASS,delta_MASS,totalparams_MASS,paramerrors_MASS] = timeseriesfit([fill_dates' A5_Centers_Sort(ins).fill_MASS'], ...
                %                 var(A5_Centers_Sort(ins).fill_MASSres),1,1);
                %             A5_Centers_Fit(ins).fill_MASS_slope=totalparams_MASS(2)*365.25;
                %             A5_Centers_Fit(ins).fill_MASS_fit = [fill_dates' fit_MASS delta_MASS];
                %             A5_Centers_Fit(ins).fill_MASS_paramerrors = paramerrors_MASS(2)*365.25; % this (365.25) is special for 'timeseriesfit' results.

                % the Slepian (MSSA) mass
                [fit_MASS,delta_MASS,totalparams_MASS,paramerrors_MASS] = timeseriesfit([fill_dates' A5_Centers_Sort(ins).fill_MASS_MSSA'], ...
                    var(A5_Centers_Sort(ins).fill_MASSres_MSSA),1,1);
                A5_Centers_Fit(ins).fill_MASS_MSSA_slope=totalparams_MASS(2)*365.25;
                A5_Centers_Fit(ins).fill_MASS_MSSA_fit = [fill_dates' fit_MASS delta_MASS];
                A5_Centers_Fit(ins).fill_MASS_MSSA_paramerrors = paramerrors_MASS(2)*365.25; % this (365.25) is special for 'timeseriesfit' results.

                % the Slepian (MSSA,accepted) mass
                % [fit_MASS,delta_MASS,totalparams_MASS,paramerrors_MASS] = timeseriesfit([fill_dates' A5_Centers_Sort(ins).fill_MASS_MSSA_both'], ...
                %     var(A5_Centers_Sort(ins).fill_MASS_MSSA_bothfail),1,1);
                [fit_MASS,delta_MASS,totalparams_MASS,paramerrors_MASS] = timeseriesfit([fill_dates' A5_Centers_Sort(ins).fill_MASS_MSSA_both'], ...
                    var(A5_Centers_Sort(ins).fill_MASS_MSSA_bothfail),1,1);
                A5_Centers_Fit(ins).fill_MASS_MSSA_both_slope=totalparams_MASS(2)*365.25;
                A5_Centers_Fit(ins).fill_MASS_MSSA_both_fit = [fill_dates' fit_MASS delta_MASS];
                A5_Centers_Fit(ins).fill_MASS_MSSA_both_paramerrors = paramerrors_MASS(2)*365.25; % this (365.25) is special for 'timeseriesfit' results.

            end

            %% EWH

            for ins=1:intit_num+1

                % the Slepian (MSSA) EWH
                [fit_EWH,delta_EWH,totalparams_EWH,paramerrors_EWH] = timeseriesfit([fill_dates' A5_Centers_Sort(ins).fill_EWH_MSSA'], ...
                    var(A5_Centers_Sort(ins).fill_EWHres_MSSA),1,1);
                A5_Centers_Fit(ins).fill_EWH_MSSA_slope=totalparams_EWH(2)*365.25;
                A5_Centers_Fit(ins).fill_EWH_MSSA_fit = [fill_dates' fit_EWH delta_EWH];
                A5_Centers_Fit(ins).fill_EWH_MSSA_paramerrors = paramerrors_EWH(2)*365.25; % this (365.25) is special for 'timeseriesfit' results.

                % the Slepian (MSSA,accepted) EWH
                [fit_EWH,delta_EWH,totalparams_EWH,paramerrors_EWH] = timeseriesfit([fill_dates' A5_Centers_Sort(ins).fill_EWH_MSSA_both'], ...
                    var(A5_Centers_Sort(ins).fill_EWH_MSSA_bothfail),1,1);
                A5_Centers_Fit(ins).fill_EWH_MSSA_both_slope=totalparams_EWH(2)*365.25;
                A5_Centers_Fit(ins).fill_EWH_MSSA_both_fit = [fill_dates' fit_EWH delta_EWH];
                A5_Centers_Fit(ins).fill_EWH_MSSA_both_paramerrors = paramerrors_EWH(2)*365.25; % this (365.25) is special for 'timeseriesfit' results.

            end

            %%
            disp('Next codes are special for some areas. Please be sure you have changed them.')
            %         pause
            %% (Here is special for some cases) compare with other products

            A5_Centers_STD_Compare=A5_Centers_STD(intit_num+1);
            A5_Centers_RMS_Compare=A5_Centers_RMS(intit_num+1);
            A5_Centers_Fit_Compare=A5_Centers_Fit(intit_num+1);

            A5_Centers_Compare=[];

            % we choose the mean of the four institutions
            A5_Centers_Compare(1).name=A5_Centers_Sort(intit_num+1).name;
            A5_Centers_Compare(1).c11cmn=c11cmn_combine;
            A5_Centers_Compare(1).lon=lon;A5_Centers_Compare(1).lat=lat;

            A5_Centers_Compare(1).data_year_beg=A5_Centers_Sort(intit_num+1).data_year_beg;A5_Centers_Compare(1).process=A5_Centers_Sort(intit_num+1).process;
            A5_Centers_Compare(1).data_dates=A5_Centers_Sort(intit_num+1).data_dates;
            A5_Centers_Compare(1).fill_dates=A5_Centers_Sort(intit_num+1).fill_dates;A5_Centers_Compare(1).fill_months=fill_months;
            A5_Centers_Compare(1).use_dates=A5_Centers_Sort(intit_num+1).use_dates;A5_Centers_Compare(1).use_months=A5_Centers_Sort(intit_num+1).use_months;
            A5_Centers_Compare(1).missing_dates=A5_Centers_Sort(intit_num+1).missing_dates;A5_Centers_Compare(1).missing_months=A5_Centers_Sort(intit_num+1).missing_months;

            % eigenvalues for signal and noise
            A5_Centers_Compare(1).fill_MSSAcoffs_eva=A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_eva;
            A5_Centers_Compare(1).fill_MSSAcoffs_RCTest=A5_Centers_Sort(intit_num+1).fill_MSSAcoffs_RCTest;

            % EWH
            A5_Centers_Compare(1).fill_EWH_Slepian=A5_Centers_Sort(intit_num+1).fill_EWH_MSSA+A5_Centers_Sort(intit_num+1).fill_EWH_MSSA_resid;
            A5_Centers_Compare(1).fill_EWH_nonoisecorrection=A5_Centers_Sort(intit_num+1).fill_EWH_MSSA;
            A5_Centers_Compare(1).fill_EWH_noisecorrection=A5_Centers_Sort(intit_num+1).fill_EWH_MSSA_bothfail;
            A5_Centers_Compare(1).fill_EWH=A5_Centers_Sort(intit_num+1).fill_EWH_MSSA_both;
            A5_Centers_Compare(1).fill_EWHsig=A5_Centers_Sort(intit_num+1).fill_EWHsig_MSSA_both;
            A5_Centers_Compare(1).fill_EWHres=A5_Centers_Sort(intit_num+1).fill_EWHres_MSSA_both;

            % uptoS
            A5_Centers_Compare(1).fill_uptoS_EWH_nonoisecorrection=A5_Centers_Sort(intit_num+1).fill_uptoS_EWH_MSSA;
            A5_Centers_Compare(1).fill_uptoS_EWH_noisecorrection=A5_Centers_Sort(intit_num+1).fill_uptoS_EWH_MSSA_bothfail;
            A5_Centers_Compare(1).fill_uptoS_EWH=A5_Centers_Sort(intit_num+1).fill_uptoS_EWH_MSSA_both;
            A5_Centers_Compare(1).fill_uptoS_EWHsig=A5_Centers_Sort(intit_num+1).fill_uptoS_EWHsig_MSSA_both;
            A5_Centers_Compare(1).fill_uptoS_EWHres=A5_Centers_Sort(intit_num+1).fill_uptoS_EWHres_MSSA_both;

            A5_Centers_Compare(1).fill_Grid_EWH_Slepian=A5_Centers_Sort(intit_num+1).fill_Grid_EWH_MSSA+A5_Centers_Sort(intit_num+1).fill_Grid_EWH_MSSA_resid;
            A5_Centers_Compare(1).fill_Grid_EWH_nonoisecorrection=A5_Centers_Sort(intit_num+1).fill_Grid_EWH_MSSA;
            A5_Centers_Compare(1).fill_Grid_EWH=A5_Centers_Sort(intit_num+1).fill_Grid_EWH_MSSA_both;
            A5_Centers_Compare(1).fill_Grid_EWHsig=A5_Centers_Sort(intit_num+1).fill_Grid_EWHsig_MSSA_both;
            A5_Centers_Compare(1).fill_Grid_EWHres=A5_Centers_Sort(intit_num+1).fill_Grid_EWHres_MSSA_both;

            % MASS
            A5_Centers_Compare(1).fill_MASS=A5_Centers_Sort(intit_num+1).fill_MASS_MSSA_both;
            A5_Centers_Compare(1).fill_MASSsig=A5_Centers_Sort(intit_num+1).fill_MASSsig_MSSA_both;
            A5_Centers_Compare(1).fill_MASSres=A5_Centers_Sort(intit_num+1).fill_MASSres_MSSA_both;

            A5_Centers_Compare(1).fill_Grid_MASS=A5_Centers_Sort(intit_num+1).fill_Grid_MASS_MSSA_both;
            A5_Centers_Compare(1).fill_Grid_MASSsig=A5_Centers_Sort(intit_num+1).fill_Grid_MASSsig_MSSA_both;
            A5_Centers_Compare(1).fill_Grid_MASSres=A5_Centers_Sort(intit_num+1).fill_Grid_MASSres_MSSA_both;

            if S_ol>S
                A5_Centers_Compare(1).fill_MSL_Slepian=A5_Centers_Sort(intit_num+1).fill_MSL_MSSA+A5_Centers_Sort(intit_num+1).fill_MSL_MSSA_resid;
                A5_Centers_Compare(1).fill_MSL=A5_Centers_Sort(intit_num+1).fill_MSL_MSSA_both;
                A5_Centers_Compare(1).fill_MSLsig=A5_Centers_Sort(intit_num+1).fill_MSLsig_MSSA_both;
                A5_Centers_Compare(1).fill_MSLres=A5_Centers_Sort(intit_num+1).fill_MSLres_MSSA_both;

                A5_Centers_Compare(1).fill_Grid_MSL_Slepian=A5_Centers_Sort(intit_num+1).fill_Grid_MSL_MSSA+A5_Centers_Sort(intit_num+1).fill_Grid_MSL_MSSA_resid;
                A5_Centers_Compare(1).fill_Grid_MSL=A5_Centers_Sort(intit_num+1).fill_Grid_MSL_MSSA_both;
                A5_Centers_Compare(1).fill_Grid_MSLsig=A5_Centers_Sort(intit_num+1).fill_Grid_MSLsig_MSSA_both;
                A5_Centers_Compare(1).fill_Grid_MSLres=A5_Centers_Sort(intit_num+1).fill_Grid_MSLres_MSSA_both;

                A5_Centers_Compare(1).fill_IBdeltacoffs=A5_Centers_Sort(intit_num+1).fill_IBdeltacoffs; %IB
                A5_Centers_Compare(1).fill_GADfilldelta=A5_Centers_Sort(intit_num+1).fill_GADfilldelta; %GAD
            end

            save(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' Attach_ALL ...
                '_TP' num2str(Turning_number) '_S' num2str(TP_ss) '_N' num2str(TP_nn) '_M' num2str(M_rec) '_Sig' num2str(Noise_SigLev) '_compare.mat']),...
                'A5_Centers_Compare','A5_Centers_STD_Compare','A5_Centers_RMS_Compare','A5_Centers_Fit_Compare',...
                'c11cmn','c11cmn_combine','Lwindow','XY','MaxNumChanges_MSSA','MaxNumChanges_V',...
                'turn_V_four','turn_MSSA_four','S','S_sig','S_ol','S_shannon','BasinArea',...
                'N_gap','N_rec','M_gap','M_rec','Turning_number','Noise_SigLev',...
                'used_sta_V','used_sta_MSSA')

        end
    end


    %% the STD of selected process

    use_TP_ss=1:2:MaxNumChanges_V;
    use_TP_nn=1:2:MaxNumChanges_MSSA;

    disp_data=4;

    test_process_num=4;
    use_process_name=["Slepian","M-SSA (residuals)","M-SSA (noises)","M-SSA (selected)"];

    for ins=1

        f6=figure;

        F_gap=[.04 .03];F_marg_h=[.06 .04];F_marg_w=[.03 .08];

        ybou1=[0 20]; % SCS
        h0=tight_subplot(numel(use_TP_ss),numel(use_TP_nn),F_gap,F_marg_h,F_marg_w);

        % You'd better adjust this manually
        set (gcf,'Position',[100,50,900,900])

        % choose nexttile or subplot (based on your matlab version containing tiledlayout or not)
        % nexttile is recommended for testing; subplot for final publish (because of the title of colorbar maybe covered by subplots)
        %     t=tiledlayout('flow'); % nexttile method
        %     [a,b]=sizewind(S); % subplot method

        ff=0;

        for TP_ss=use_TP_ss

            SSF=turn_V_four(TP_ss,used_sta_V);
            % choose nexttile or subplot (based on your matlab version containing tiledlayout or not)
            %             nexttile % nexttile method
            %         subplot(a,b,SSF)  % subplot method

            for TP_nn=use_TP_nn

                load(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' ...
                    Attach_ALL '_TP' num2str(Turning_number) '_S' num2str(TP_ss) ...
                    '_N' num2str(TP_nn) '_M' num2str(M_rec)  '_Sig' num2str(Noise_SigLev) ...
                    '_compare.mat']),...
                    'A5_Centers_STD_Compare');

                ff=ff+1;

                axes(h0(ff));

                switch disp_data
                    case 1
                        r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_uptoS_STD(SSF,:,:));fig_name=['fill_Grid_EWH_MSSA_uptoS_STD'];
                    case 2
                        r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_resid_uptoS_STD(SSF,:,:));fig_name=['fill_Grid_EWH_MSSA_resid_uptoS_STD'];
                    case 3
                        r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_bothfail_uptoS_STD(SSF,:,:));fig_name=['fill_Grid_EWH_MSSA_bothfail_uptoS_STD'];
                    case 4
                        r_res=squeeze(A5_Centers_STD_Compare(ins).fill_Grid_EWH_MSSA_both_uptoS_STD(SSF,:,:));fig_name=['fill_Grid_EWH_MSSA_both_uptoS_STD'];
                end

                %
                m_proj('mercator','long',[c11cmn_combine(1)-1.5, c11cmn_combine(3)+0.5],'lat',[c11cmn_combine(4)-0.5, c11cmn_combine(2)+1.5]);
                m_pcolor(lonlon,latlat,r_res);
                % m_coast('patch',[1 .85 .7]);
                m_coast
                m_grid('box','fancy','tickdir','in');
                m_line(XY(:,1),XY(:,2),'color','k','linewidth',1);
                % xlabel(['CC(' num2str(SSF) ')'])

                m_text(c11cmn_combine(1)+1,c11cmn_combine(2)-1,['S=S(' num2str(TP_ss) ') N=N(' num2str(TP_nn) ')'],'color','red','fontsize',12);
                %             title(['N = ' char(use_N_name(SSF_n))],'color','red','fontsize',12)

                %             if SSF_n==1
                %                 title(use_process_name(nn),'fontsize',14)
                %             end
                %         caxis([0 10]);
                caxis(ybou1);

                set(gca,'FontName','Times New Roman','FontSize',font_Size-2,'TickLength',[0.004,0.035]);

%                 writegeofiles(lon,lat,r_res,fullfile(txt_path_ALL,[fig_name '_TP_S' num2str(TP_ss) '_N' num2str(TP_nn) '.txt']),[],4)

            end

            %%% nexttile method
            %         t.TileSpacing='Compact';
            %         t.Padding='Compact';

            colormap(jet)
            %             colormap(flipud (jet))
            %         colormap(mymap ('Greys'))

        end

        cb = colorbar;
        cb.Position = [0.94 0.1 0.02,0.8];
        title(cb,'STD','fontsize',14);
        %%%

        whratio=abs((c11cmn_combine(3)-c11cmn_combine(1))/c11cmn_combine(2)-c11cmn_combine(4));

        tif_name6=[fig_name,'_Spatial_STD_NSchange.tif'];
        print(f6,'-dtiff','-r300',fullfile(char(fig_path_ALL),tif_name6));


    end

end