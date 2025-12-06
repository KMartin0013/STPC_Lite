% This is the main code to calcualte the spherical Slepian Functions,
% enough for advance users and also friendly to beginners.
% Code version: 1 (2nd round RSE revision).
% Code Note:
% Main function:
%   (1) calculate the total mass change (in Gt) for the given basin, and
%   the corresonding spherical Slepian Functions for further process.
% Main features:
%   (1) fill the gap within and among GRACE and GRACE-FO datasets by
%   eight-coefficient fits (Gauer et al., 2022). This fitted values will
%   serve as the raw data for iterative M-SSA gap fitting procedure.
%   (2) consider the case of artificial gap for further sensitivity
%   analysis.
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
load('Slepian_Information');

%% Main code
%% Projects GRACE into the requested Slepian basis

FIG_Attach=sprintf('%s_%s_%s_%s_%s_%s',...
    Dataproduct{1},Dataproduct{2}(3:4),Area,num2str(Lwindow),buffer_str,num2str(Radius));
fnpl_IB=sprintf('%s/%s_%s_IB_%s_%s_%s.mat',fullfile(getenv('IFILES'),'GRACE'),...
    Dataproduct{1},Dataproduct{2},num2str(Lwindow),'SD','GSHHC');

load([getenv('IFILES'), '/EARTHMODELS/CONSTANTS/Earth.mat']);
%     radius=Earth.Radius;

Earth_radius=fralmanac('a_EGM96','Earth'); % be consistent with code
%     Earth_radius=6370000; % be consistent with code

%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0=figure;
XY_buf=eval(sprintf('%s(%i,%f)',TH_ori,0,buffer_deg));
BasinArea_buffer=spharea(XY_buf)*4*pi*Earth_radius^2;
geoshow(XY_buf(:,2), XY_buf(:,1), 'DisplayType', 'polygon', ...
    'FaceColor', 'yellow')

hold on

XY_ori=eval(sprintf('%s(%i,%f)',TH_ori,0,0));
BasinArea_ori=spharea(XY_ori)*4*pi*Earth_radius^2;

% this is too empirical, so the buffer contour and study area my be
% separated with 360°
if XY_ori(1,1)>180
    plot(XY_ori(:,1)-360, XY_ori(:,2),'--b','linewidth',3)
else
    plot(XY_ori(:,1), XY_ori(:,2),'--b','linewidth',3)
end

% plot(XY_ori(:,1),XY_ori(:,2),'b','linewidth',3);

title('Study area');

tif_name0=['Fig0_',FIG_Attach,'.tif'];
print(f0,'-dtiff','-r125',fullfile(getenv('IFILES'),['Figure_' Code_Version],tif_name0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mass density
[data_slepcoffs,calerrors,data_dates,TH_buffer,G,CC,V]=grace2slept_m(Dataproduct,...
    TH_ori,buffer_deg,Lwindow,phi,theta,omega,[],'SD',1);

% choose the final basin for quantitative analysis
%   This is for the whole land without counting for the buffer zone becuase
%   the land signal generally leaking to nearby ocean while suitable
%   selection of buffer zone 'recover' this leakage. However, in specific
%   cases (such as greenland), the leakage is not only due to the filter
%   process, but also the truncation of GRACE (up to 60/96 degree). Therefore
%   some grids around the the land should be accounted.
%
%   For ocean, however,
%   buffer zone should be considered to avoid the strong leakage by land.
if strcmp(Dataproduct{4},'land') && ~strcmp(TH_ori,'greenland')
    TH=TH_ori;XY=XY_ori;
    BasinArea=BasinArea_ori;
    %     TH=TH_buffer;XY=XY_buf;
    %     BasinArea=BasinArea_buffer;
else
    TH=TH_buffer;XY=XY_buf;
    BasinArea=BasinArea_buffer;
end

save(fullfile(ddir1,['XY_',FIG_Attach,'.mat']),'XY_buf','XY_ori','XY')

S_shannon=round((Lwindow+1)^2*spharea(XY_buf));
switch S_choice
    case 0
        % contain Shannon number
        defval('S',round((Lwindow+1)^2*spharea(XY_buf)));
        S_sig=S;
    case 1
        % contain signal > 0.1
        S=sum(V>0.1);
        S_sig=S;

    case 2
        % contain also additional smooth
        S=sum(V>0.1);
        N1=round((Lwindow+1)^2*spharea(XY_buf));
        N2=sum(V>0.1);
        S_sig=[N1,N2];

    case 3
        % smooth all SSF
        defval('S',round((Lwindow+1)^2*spharea(XY_buf)));
        N1=0;
        N2=S;
        S_sig=[N1,N2];
    case 4
        % contain signal > 0.3
        S=sum(V>0.3);
        S_sig=S;
    case 5
        % contain signal > 0.01
        S=sum(V>0.01);
        S_sig=S;
    case 6
        % S=50
        S=Max_S;
        S_sig=S;

    otherwise
        % default use Shannon number
        disp('Given no preference on S, the shannon number is selected.')
end

% defval('S',50);
if length(S_sig)>1
    FIG_Attach=[FIG_Attach '_S' num2str(N1) 't' num2str(N2)]
else
    FIG_Attach=[FIG_Attach '_S' num2str(S)]
end

%% Additional code for finding leakaging date
% this is a automatic process and normally free of manual adjustment

InsYear_str=char(Dataproduct(1));
% filled time epochs for 2005-2023 without gaps
fill_nmonths=(str2num(InsYear_str(end-1:end))-str2num(InsYear_str(end-3:end-2))+1)*12;

data_year_beg=2000+str2num(InsYear_str(end-3:end-2));
data_year_end=2000+str2num(InsYear_str(end-1:end));

% raw data month (may include repeat months)
data_months=(str2num(datestr(data_dates,'yyyy'))-data_year_beg)*12+ ...
    str2num(datestr(data_dates,'mm'));
data_nmonths=numel(data_months);

% Transfer from data_months to use_months (without repeat months)
%   Because in some GRACE products there may be more than one value for a
%   month, this procedure is designed to check for any possible duplicate
%   months. Therefore, the 'data_mon' is not always the same to the
%   'use_months' afterwards.
[use_months,data2use_monIndex] = unique(data_months,'stable'); % IndexTimeA=
use_nmonths = length(use_months);use_dates=data_dates(data2use_monIndex);

RepeatItem = 0; % The first number, 0, is invalid for indexing purposes
for ii = 2:use_nmonths % find duplicates
    if data2use_monIndex(ii) - data2use_monIndex(ii-1) > 1
        NowRepeatItem = [0,data2use_monIndex(ii-1):(data2use_monIndex(ii)-1)];
        RepeatItem = [RepeatItem,NowRepeatItem];
    end
end
if numel(RepeatItem)>1
    disp('Duplicate months appear in')
    disp(RepeatItem(2:end));
end

%% Slepian residuals and integration
% 3. Next run SLEPT2RESID to fit a choice of functions (e.g. lines,
% quadratics, etc.) to the Slepian coefficients. Note: if you want to
% remove a model of GIA then you should do that before this step, using
% a function like CORRECT4GIA.

%     fitwhat=[3 365.25]; %quadratic periodic residual and an annual harmonic term
fitwhat=[3 365.25 365.25/2]; %linear trend and an annual and a semiannual harmonic term

% 4. If you want to examine the total mass trend in a region, this
% information is calculated and returned as a result from SLEPT2RESID.
% To summarize, each Slepian coefficient up to the Shannon truncation is
% multiplied by the integral (found using INTEGRATEBASIS) of the corresponding
% function over the region. This total data is then fit with TIMESERIESFIT
% which allows fitting with a custom data variance.

% transfer use_mon to EST_mon for further gap interpolation
fill_months=[1:fill_nmonths]';

% find the missing time epoch (gap) in and among GRACE and GRACE-FO
missing_months=setdiff(fill_months,use_months)';
Slept_months=[data_months;missing_months'];
if any(missing_months)

    %fill the time epoch by linear interpolate (abandon)
    %     thedates_int = interp1(given_mon,thedates,leakage);
    %     ESTthedates = interp1(given_mon,thedates,EST_mon)';

    %fill the time epoch by the mid-month date of each leakaging month
    for i=1:numel(missing_months)
        missing_monthstart = datenum([data_year_beg+floor((missing_months(i)-1)/12) ...
            mod(missing_months(i)-1,12)+1 1]);
        missing_monthend = datenum([data_year_beg+floor((missing_months(i))/12) ...
            mod(missing_months(i),12)+1 1])-1;
        missing_monthmid = (missing_monthstart+missing_monthend)/2;
        missing_dates(i) = missing_monthmid;
    end
    Slept_dates=[data_dates missing_dates];

    fill_dates=zeros(1,numel(fill_months));
    fill_dates(use_months)=use_dates;
    fill_dates(missing_months)=missing_dates;
else

    Slept_dates=data_dates;
    fill_dates=data_dates;
end

% finally, calculate the slepian functions, fit the slepian coefficients,
% and get the signal and residual for use_months and also interpolated
% signal for leakage_months. And also, get the total change over the study
% area.
[data_signalcoffs,data_residcoffs,data_ftests,data_extravalues,data_MASS,...
    data_alphavarall,data_MASSparams,data_MASSparamerrors,...
    data_MASSfit,functionintegrals,alphavar]=slept2resid_m(data_slepcoffs,...
    Slept_dates,fitwhat,[],[],CC,TH,S_sig,Radius);

% try my best to be consistent with previous definitions that each month
% should have one and only have one value.
% slepian coefficients (unit: kg/m^2)
use_signalcoffs=data_signalcoffs(data2use_monIndex,:);
use_residcoffs=data_residcoffs(data2use_monIndex,:);
% total mass change and fit (unit: Gt)
use_MASS=data_MASS(data2use_monIndex);
use_MASSfit=data_MASSfit(data2use_monIndex,:);


% this codes are for special sensitivity analysis by summing a artificial
% gap in the sixth year of your study area. Notice that you had better
% garantee there is no missing month in this period.
if exist('artificial_months','var')

    %     artificial_months=5*12+1:6*12;

    Art_data_months=data_months;Art_data_dates=data_dates;
    Art_data_slepcoffs=data_slepcoffs;
    [ia,ib] = ismember(Art_data_months,artificial_months);
    Art_data_months(ia)=[];Art_data_dates(ia)=[];
    Art_data_slepcoffs(ia,:)=[];

    Art_use_months=use_months;Art_use_dates=use_dates;
    [ia,ib] = ismember(Art_use_months,artificial_months);
    Art_use_months(ia)=[];Art_use_dates(ia)=[];
    Art_use_nmonths = length(Art_use_months);

    % raw data month (may include repeat months)
    Art_data_months=(str2num(datestr(Art_data_dates,'yyyy'))-data_year_beg)*12+ ...
        str2num(datestr(Art_data_dates,'mm'));
    Art_data_nmonths=numel(Art_data_months);

    % find repeat month (if repeat, use the first values)
    % transfer from data_mon to use_mon (without repeat months)
    [Art_use_months,Art_data2use_monIndex] = unique(Art_data_months,'stable'); % IndexTimeA返回索引

    for i=1:numel(artificial_months)
        artificial_monthstart = datenum([data_year_beg+floor((artificial_months(i)-1)/12) ...
            mod(artificial_months(i)-1,12)+1 1]);
        artificial_monthend = datenum([data_year_beg+floor((artificial_months(i))/12) ...
            mod(artificial_months(i),12)+1 1])-1;
        artificial_monthmid = (artificial_monthstart+artificial_monthend)/2;
        artificial_dates(i) = artificial_monthmid;
    end
    Art_Slept_dates=[Art_data_dates missing_dates artificial_dates];


    [Art_data_signalcoffs,Art_data_residcoffs,Art_data_ftests,Art_data_extravalues]=slept2resid_m(Art_data_slepcoffs,...
        Art_Slept_dates,fitwhat,[],[],CC,TH,S_sig,Radius);

    % try my best to be consistent with previous definitions that each month
    % should have one and only have one value.
    Art_ESTsignal=Art_data_signalcoffs(Art_data2use_monIndex,:);
    Art_ESTresid=Art_data_residcoffs(Art_data2use_monIndex,:);

end
%% Temporal analysis on the Slepian coefficients
%     given_nmonths=length(thedates); % given_nmonths ~ use_mon_num

fill_nmonths=length(fill_dates);

[Cab] = slepresid2cov(use_residcoffs);

% Make the coefficients with reference to some mean
% If they already are, then this won't matter
use_slepcoffs=data_slepcoffs(data2use_monIndex,:); % without repeat months
% the mean is consistent with the inner code in slept2resid_m (include repeat months)
use_sleptdelta = use_slepcoffs(1:use_nmonths,:) - repmat(mean(data_slepcoffs(1:data_nmonths,:),1),use_nmonths,1);

[use_MASS_CC,~,~,~,use_MASSfit_CC,~]=integral_fit(S,use_dates,functionintegrals,use_sleptdelta,Cab);

if any(missing_months)
    % complete estiamted signal with filled leakage values
    fill_signalcoffs_com=zeros(length(fill_dates),size(use_signalcoffs,2));
    fill_residcoffs_com=zeros(length(fill_dates),size(use_residcoffs,2));
    for i=1:size(use_signalcoffs,1)
        fill_signalcoffs_com(use_months(i),:)=use_signalcoffs(i,:);
        fill_residcoffs_com(use_months(i),:)=use_residcoffs(i,:);
    end

    %     extravalues_resid=zeros(length(missing_months),size(ESTresid,2));
    %     % old version for gap interpolate
    %     for j=1:size(ESTresid_com,2)
    %         ESTresid_missing = interp1(use_dates,ESTresid(:,j)',missing_dates,'spline');
    %         extravalues_resid(:,j)=ESTresid_missing';
    %     end

    for i=1:size(missing_months,2)
        fill_signalcoffs_com(missing_months(i),:)=data_extravalues(i,:);
        %             ESTresid_com(leakage(i),:)=extravalues_resid(i,:);
        fill_residcoffs_com(missing_months(i),:)=0;
    end

    fill_signaldelta=fill_signalcoffs_com(1:fill_nmonths,:) - repmat(mean(use_signalcoffs(1:use_nmonths,:),1),fill_nmonths,1);
    fill_residdelta=fill_residcoffs_com(1:fill_nmonths,:) - repmat(mean(use_residcoffs(1:use_nmonths,:),1),fill_nmonths,1);

    if any(artificial_months)
        % complete estiamted signal with filled leakage values
        Art_fill_signalcoffs_com=zeros(length(fill_dates),size(Art_ESTsignal,2));
        Art_fill_residcoffs_com=zeros(length(fill_dates),size(Art_ESTresid,2));
        for i=1:size(Art_ESTsignal,1)
            Art_fill_signalcoffs_com(Art_use_months(i),:)=Art_ESTsignal(i,:);
            Art_fill_residcoffs_com(Art_use_months(i),:)=Art_ESTresid(i,:);
        end

        %         extravalues_missing_resid=zeros(length(missing_months),size(Art_ESTresid,2));
        %         extravalues_art_resid=zeros(length(artificial_months),size(Art_ESTresid,2));
        %         % old version for gap interpolate
        %         for j=1:size(Art_ESTresid_com,2)
        %             Art_ESTresid_missing = interp1(Art_use_dates,Art_ESTresid(:,j)',missing_dates,'spline');
        %             extravalues_missing_resid(:,j)=Art_ESTresid_missing';
        %
        %             Art_ESTresid_art = interp1(Art_use_dates,Art_ESTresid(:,j)',artificial_dates,'spline');
        %             extravalues_art_resid(:,j)=Art_ESTresid_art';
        %         end

        for i=1:size(missing_months,2)
            Art_fill_signalcoffs_com(missing_months(i),:)=Art_data_extravalues(i,:);
            %             ESTresid_com(leakage(i),:)=extravalues_resid(i,:);
            Art_fill_residcoffs_com(missing_months(i),:)=0;
        end

        for i=1:size(artificial_months,2)
            Art_fill_signalcoffs_com(artificial_months(i),:)=Art_data_extravalues(numel(missing_months)+i,:);
            %             ESTresid_com(leakage(i),:)=extravalues_resid(i,:);
            Art_fill_residcoffs_com(artificial_months(i),:)=0;
        end

        Art_fill_signaldelta=Art_fill_signalcoffs_com(1:fill_nmonths,:) - repmat(mean(Art_ESTsignal(1:Art_use_nmonths,:),1),fill_nmonths,1);
        Art_fill_residdelta=Art_fill_residcoffs_com(1:fill_nmonths,:) - repmat(mean(Art_ESTresid(1:Art_use_nmonths,:),1),fill_nmonths,1);
    end

else
    fill_signaldelta=use_signalcoffs(1:use_nmonths,:) - repmat(mean(use_signalcoffs(1:use_nmonths,:),1),use_nmonths,1);
    fill_residdelta=use_residcoffs(1:use_nmonths,:) - repmat(mean(use_residcoffs(1:use_nmonths,:),1),use_nmonths,1);

    if any(artificial_months)
        % complete estiamted signal with filled leakage values
        Art_fill_signalcoffs_com=zeros(length(fill_dates),size(Art_ESTsignal,2));
        Art_fill_residcoffs_com=zeros(length(fill_dates),size(Art_ESTresid,2));
        for i=1:size(Art_ESTsignal,1)
            Art_fill_signalcoffs_com(Art_use_months(i),:)=Art_ESTsignal(i,:);
            Art_fill_residcoffs_com(Art_use_months(i),:)=Art_ESTresid(i,:);
        end

        %         extravalues_art_resid=zeros(length(artificial_months),size(Art_ESTresid,2));
        %         % old version for gap interpolate
        %         for j=1:size(Art_fill_residcoffs_com,2)
        %             Art_ESTresid_art = interp1(Art_use_dates,Art_ESTresid(:,j)',artificial_dates,'spline');
        %             extravalues_art_resid(:,j)=Art_ESTresid_art';
        %         end

        for i=1:size(artificial_months,2)
            Art_fill_signalcoffs_com(artificial_months(i),:)=Art_data_extravalues(i,:);
            %             ESTresid_com(leakage(i),:)=extravalues_resid(i,:);
            Art_fill_residcoffs_com(artificial_months(i),:)=0;
        end

        Art_fill_signaldelta=Art_fill_signalcoffs_com(1:fill_nmonths,:) - repmat(mean(Art_ESTsignal(1:Art_use_nmonths,:),1),fill_nmonths,1);
        Art_fill_residdelta=Art_fill_residcoffs_com(1:fill_nmonths,:) - repmat(mean(Art_ESTresid(1:Art_use_nmonths,:),1),fill_nmonths,1);
    end

end

% this is the filled spherical Slepian functions
fill_deltacoffs=fill_signaldelta+fill_residdelta; % unit kg/m^2

if any(artificial_months)
    Art_fill_delta=Art_fill_signaldelta+Art_fill_residdelta; % unit kg/m^2
end

% calculate the mass change for each SSF
% notion that here Cab could not be used actually. So we do not use the
% error calculation here.
[fill_MASSsig_CC,~,~,~,fill_MASSsigfit_CC,~]=integral_fit(S,fill_dates,functionintegrals,fill_signaldelta,Cab);

% notion that here Cab could not be used actually. So we do not use the
% error calculation here.
[fill_MASSres_CC,~,~,~,fill_MASSresfit_CC,~]=integral_fit(S,fill_dates,functionintegrals,fill_residdelta,Cab);

fill_MASSsig=sum(fill_MASSsig_CC);
fill_MASSres=sum(fill_MASSres_CC);
fill_MASS=sum(fill_MASSsig_CC)+sum(fill_MASSres_CC);

% figure
% plot(use_months,use_MASS);
% hold on
% plot(fill_months,fill_MASS);
% legend('UseMass','FillMass')
%% transfer to Equivelent water height (spatial)
% this Equavalent water height could be terrestrial water storage or Ocean
% bottom pressure depanding on your study area in land or ocean
% unit: cm

[r_example,lon,lat,Plm,degres]=plm2xyz(CC{1},1,c11cmn,Lwindow);

if lon<0
    lon=lon+360;
end

% % area-weighted can be used to more precisely calculate the regional
% % basin-averaged variables
[lonlon,latlat]=meshgrid(lon,lat);
[in, on]=check_polygon_in(XY,lonlon,latlat);
for i=1:length(lat)
    for j=1:length(lon)
        c11cmn_area(i,j)=areaquad(lat(i)-0.5,lon(j)-0.5,lat(i)+0.5,lon(j)+0.5);
    end
end

if length(S_sig)>1

    for j=1:N1
        [r_record(j,:,:),lon,lat,Plm,degres]=plm2xyz(CC{j},1,c11cmn,Lwindow);
    end
    for j=N1+1:N2
        CC_smo=plm_Gausssmooth(CC{j},Radius);
        [r_record(j,:,:),lon,lat,Plm,degres]=plm2xyz(CC_smo,1,c11cmn,Lwindow);
    end
else

    for j=1:S
        [r_record(j,:,:),lon,lat,Plm,degres]=plm2xyz(CC{j},1,c11cmn,Lwindow);
    end
end

fill_awMASS=zeros(1,numel(fill_dates));
fill_Grid_EWHsig=zeros(fill_nmonths,size(r_example,1),size(r_example,2));
fill_Grid_EWHres=zeros(fill_nmonths,size(r_example,1),size(r_example,2));
for i=1:fill_nmonths
    sp_ewh_signal=0;sp_ewh_residual=0;
    for j=1:S
        %             [r,lon,lat,Plm,degres]=plm2xyz(CC{j},1,c11cmn,Lwindow);
        r=squeeze(r_record(j,:,:));

        sp_ewh_signal=sp_ewh_signal+r*fill_signaldelta(i,j)' ... % (you can manually adjust this period)
            /1000*100; % change from kg/m^2 to cm for pure water (/ 1000 kg/m*3 * 100)
        sp_ewh_residual=sp_ewh_residual+r*fill_residdelta(i,j)' ... % (you can manually adjust this period)
            /1000*100; % change from kg/m^2 to cm for pure water (/ 1000 kg/m*3 * 100)
    end
    % substract IB
    fill_Grid_EWHsig(i,:,:)=sp_ewh_signal;
    fill_Grid_EWHres(i,:,:)=sp_ewh_residual;

    fill_awMASS(i)=sum( (sp_ewh_signal(in)+sp_ewh_residual(in))/100.*c11cmn_area(in) )*4*pi*6370000^2*10^3/10^3/10^9; % Gt
end

fill_Grid_EWH=fill_Grid_EWHsig+fill_Grid_EWHres;

figure
plot(1:fill_nmonths,fill_MASS);
hold on
plot(1:fill_nmonths,fill_awMASS);

legend('Code','Area-weighted');


% signal + residual (raw data)
use_Grid_EWH=zeros(use_nmonths,size(r_example,1),size(r_example,2));
for i=1:use_nmonths
    sp_ewh=0;
    for j=1:S
        %             [r,lon,lat,Plm,degres]=plm2xyz(CC{j},1,c11cmn,Lwindow);
        r=squeeze(r_record(j,:,:));

        sp_ewh=sp_ewh+r*use_sleptdelta(i,j)' ... % (you can manually adjust this period)
            /1000*100; % change from kg/m^2 to cm (/ 1028 kg/m*3 * 100)
    end
    use_Grid_EWH(i,:,:)=sp_ewh;
end

save(fullfile(ddir1,['MainData_' FIG_Attach]),"Dataproduct",...
    "use_sleptdelta","use_residcoffs","use_signalcoffs","use_MASS_CC", ...
    "use_MASS","use_MASSfit","use_MASSfit_CC", ...
    "functionintegrals","Cab",...
    "fill_signaldelta","fill_residdelta","fill_deltacoffs",...
    "fill_MASS","fill_MASSsig","fill_MASSres",...
    "fill_MASSsig_CC","fill_MASSsigfit_CC",...
    "fill_MASSres_CC","fill_MASSresfit_CC",...
    "fill_Grid_EWH","fill_Grid_EWHsig","fill_Grid_EWHres",...
    "data_MASSparams","data_MASSparamerrors","data_alphavarall",...
    "data_year_beg","use_dates","fill_dates","Slept_dates","missing_dates", ...
    "use_months","fill_months","Slept_months","missing_months",...
    "data2use_monIndex","data_dates","data_months",...
    "CC","V","c11cmn","S","S_sig","S_shannon","Radius", ...
    "TH","BasinArea","XY","buffer_str","Area","Lwindow")


% Artificial
if any(artificial_months)
    save(fullfile(ddir1,['MainData_' FIG_Attach]),'-append',"Art_fill_delta",...
        "artificial_months")
end

%% if you want to study in ocean, you need to correct IB this is for regional SL

if strcmp(Dataproduct{4},'ocean')

    %     IB=D_choose_SL(13);
    load(fnpl_IB)
    save(fnpl_IB,'IB_C','data_dates','use_months')
    %     IB=smooth(IB_C(:,2)/10,9); % no need to smooth here
    data_IB=smooth(IB_C(:,2),1); % unit kg/m^2 or mmH2O
    % recommend to smooth it by 90 days (i.e., about 3 months)

    % consider IB in error propagation

    % only use the fit function in this code to get the residual of IB
    [data_IBsignal,data_IBresid,data_IB_ftests,data_IB_extravalues]=slept2resid_m(data_IB,...
        Slept_dates,fitwhat); % mm (values equivalent to kg/m^2)

    use_IBsignal=data_IBsignal(data2use_monIndex,:);
    use_IBresid=data_IBresid(data2use_monIndex,:);

    % here we calculate the pseudo ocean mass correcting the IB. Notice
    % that this correction is only used for calculate the mass sea level
    % changes, DO NOT represent the actual correction on the total mass of
    % the ocean. If you want to represent the total mass loss of the ocean,
    % you should still use fill_MASS instead.

    % substitude the S+1 slepian function as IB, because we only add up slepian
    % functions up to S, the S+1 was used to correct the IB mannualy
    use_sleptdelta_includeIB=use_sleptdelta;
    use_sleptdelta_includeIB(:,S+1)=data_IB(data2use_monIndex);

    % error propagation
    use_residcoffs_includeIB=use_residcoffs;
    use_residcoffs_includeIB(:,S+1)=use_IBresid;
    [Cab_includeIB] = slepresid2cov(use_residcoffs_includeIB);
    functionintegrals_includeIB=[functionintegrals -1/10^3*BasinArea*10^3/10^3/10^9]; % subtract IB (change unit from mm to Gt)
    alphavarall_includeIB=functionintegrals_includeIB*Cab_includeIB(1:S+1,1:S+1)*functionintegrals_includeIB';

    % After correcting IB, this should be regional mass sea level
    [use_MASScIB_CC]=integral_fit(S+1,use_dates,functionintegrals_includeIB,use_sleptdelta_includeIB,Cab_includeIB);

    use_MASScIB=sum(use_MASScIB_CC);
    [MASScIB_fit,MASScIB_delta,MASScIB_totalparams,MASScIB_paramerrors] = timeseriesfit([use_dates' use_MASScIB'],...
        alphavarall_includeIB,1,1);
    % Make a matrix for the line, and 95% confidence in the fit
    use_MASScIBfit = [use_dates' MASScIB_fit MASScIB_delta];
    use_MASScIBparamerrors = MASScIB_paramerrors*365.25;

    if any(missing_months)
        % complete estiamted signal with filled leakage values
        fill_IBsignal_com=zeros(length(fill_dates),size(use_IBsignal,2));
        fill_IBresid_com=zeros(length(fill_dates),size(use_IBresid,2));
        for i=1:size(use_signalcoffs,1)
            fill_IBsignal_com(use_months(i),:)=use_IBsignal(i,:);
            fill_IBresid_com(use_months(i),:)=use_IBresid(i,:);
        end

        %     extravalues_IBresid=zeros(length(missing_months),size(IB_resid,2));
        %     % old version for gap interpolate
        %     IBresid_missing = interp1(use_dates,IB_resid(:,1)',missing_dates,'spline');
        %     extravalues_IBresid(:,1)=IBresid_missing';

        for i=1:size(missing_months,2)
            fill_IBsignal_com(missing_months(i),:)=data_IB_extravalues(i,:);
            %             IBresid_com(leakage(i),:)=extravalues_IBresid(i,:);
            fill_IBresid_com(missing_months(i),:)=0;
        end
        fill_IBsignaldelta=fill_IBsignal_com(1:fill_nmonths,:) - repmat(mean(use_IBsignal(1:use_nmonths,:),1),fill_nmonths,1);
        fill_IBresiddelta=fill_IBresid_com(1:fill_nmonths,:) - repmat(mean(use_IBresid(1:use_nmonths,:),1),fill_nmonths,1);

    else
        fill_IBsignaldelta=use_IBsignal(1:use_nmonths,:) - repmat(mean(use_IBsignal(1:use_nmonths,:),1),use_nmonths,1);
        fill_IBresiddelta=use_IBresid(1:use_nmonths,:) - repmat(mean(use_IBresid(1:use_nmonths,:),1),use_nmonths,1);

    end

    % this is the filled spherical Slepian functions
    fill_IBdeltacoffs=fill_IBsignaldelta+fill_IBresiddelta; % unit kg/m^2

    %
    %     fill_MASSsig_IB=fill_IBsignaldelta'/10^3*BasinArea*10^3/10^3/10^9; % change unit of IB from mmH2O to Gt
    %     fill_MASSres_IB=fill_IBresiddelta'/10^3*BasinArea*10^3/10^3/10^9; % change unit of IB from mmH2O to Gt
    %     fill_MASS_IB=fill_MASSsig_IB+fill_MASSres_IB;

    % %     MSL_ESTsigtotal=sum(ESTsigtotal_CC)-IBsignal';
    %     fill_MASSsig_Ocean=fill_MASSsig-fill_IBsignaldelta'/10^3*BasinArea*10^3/10^3/10^9; % change unit of IB from mmH2O to Gt

    %     fill_MASSres_Ocean=fill_MASSres-fill_IBresiddelta'/10^3*BasinArea*10^3/10^3/10^9; % change unit of IB from mmH2O to Gt

    %     fill_MASS_Ocean=fill_MASSsig_Ocean+fill_MASSres_Ocean;

    %% transfer to Equivelent water height (spatial)
    % this Equavalent water height could be terrestrial water storage or Ocean
    % bottom pressure depanding on your study area in land or ocean

    fill_Grid_MSLsig=zeros(fill_nmonths,size(r_example,1),size(r_example,2));
    fill_Grid_MSLres=zeros(fill_nmonths,size(r_example,1),size(r_example,2));
    %     fill_OBPsig=zeros(fill_nmonths,size(r_example,1),size(r_example,2));
    %     fill_OBPres=zeros(fill_nmonths,size(r_example,1),size(r_example,2));
    for i=1:fill_nmonths
        sp_obp_signal=0;sp_obp_residual=0;
        for j=1:S
            %             [r,lon,lat,Plm,degres]=plm2xyz(CC{j},1,c11cmn,Lwindow);
            r=squeeze(r_record(j,:,:));

            sp_obp_signal=sp_obp_signal+r*fill_signaldelta(i,j)' ... % (you can manually adjust this period)
                /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
            sp_obp_residual=sp_obp_residual+r*fill_residdelta(i,j)' ... % (you can manually adjust this period)
                /1000*100; % change from kg/m^2 to cm (/ 1000 kg/m*3 * 100)
        end
        % substract IB
        fill_Grid_MSLsig(i,:,:)=sp_obp_signal*1000/1028-fill_IBsignaldelta(i)'*1000/1028/10; % change unit from mmH2O to cm(salty H2O)
        fill_Grid_MSLres(i,:,:)=sp_obp_residual*1000/1028-fill_IBresiddelta(i)'*1000/1028/10; % change unit from mmH2O to cm(salty H2O)
        %         fill_OBPsig(i,:,:)=sp_obp_signal;
        %         fill_OBPres(i,:,:)=sp_obp_residual;
    end

    fill_Grid_MSL=fill_Grid_MSLsig+fill_Grid_MSLres;
    %     fill_OBPreconst=fill_OBPsig+fill_OBPres;

    % signal + residual (raw data)
    use_Grid_MSL=zeros(use_nmonths,size(r_example,1),size(r_example,2));
    %     use_OBP=zeros(use_nmonths,size(r_example,1),size(r_example,2));
    for i=1:use_nmonths
        sp_obp=0;
        for j=1:S
            %             [r,lon,lat,Plm,degres]=plm2xyz(CC{j},1,c11cmn,Lwindow);
            r=squeeze(r_record(j,:,:));

            sp_obp=sp_obp+r*use_sleptdelta(i,j)' ... % (you can manually adjust this period)
                /1000*100; % change from kg/m^2 to cm (/ 1028 kg/m*3 * 100)
        end
        % add IB
        use_Grid_MSL(i,:,:)=sp_obp*1000/1028-data_IB(data2use_monIndex(i))'*1000/1028/10; % change unit from mmH2O to cm(salty H2O)
        %         use_OBP(i,:,:)=sp_obp;
    end

    save(fullfile(ddir1,['MainData_' FIG_Attach]),'-append',...
        "fill_Grid_MSLsig","fill_Grid_MSLres","fill_Grid_MSL",...
        "fill_IBsignaldelta","fill_IBdeltacoffs",...
        "fill_IBresiddelta",...
        "use_Grid_MSL","data_IB")

end

disp('Done. Save files are in /Results.')
