function extract_GAD_degree1(Pcenter,Rlevel)

defval('Pcenter','CSR')
defval('Rlevel','RL06')

%m in case different time span of data is used
if length(Pcenter)>=3
    Pcenters=Pcenter(1:3);
end

% Where the original data files are kept
% defval('ddir1',fullfile('D:','GRACE_DATA',[Pcenter Rlevel]));
defval('ddir1',fullfile(getenv('IFILES'),'GRACE','Originals',[Rlevel '0324'])); 
% Where you would like to save the new file
defval('ddir2',fullfile(getenv('IFILES'),'GRACE','Degree1'));

if ~exist(ddir1, 'dir')
    error('Please check you route for the GAD products of GRACE. You had better change the ddir1 accordingly.')
end

ouput_filename=fullfile(ddir2,['GAD_Degree1_GEOC_' Pcenters '_' Rlevel '.txt']);

if strcmp(Pcenters,'GFZ')
    if strcmp(Rlevel,'RL04')
        % Find the coefficient files
        datanames=ls2cell(fullfile(ddir1,'GSM*G---_0004'));
        % Know a priori what the bandwidth of the coefficients is
        Ldata=120;
    elseif strcmp(Rlevel,'RL05')
        % Find the coefficient files
        datanames=ls2cell(fullfile(ddir1,'GAD*G---_0005'));
        % Know a priori what the bandwidth of the coefficients is
        Ldata=100;
    elseif strcmp(Rlevel,'RL06')
        % Find the coefficient files
        datanames=ls2cell(fullfile(ddir1,'GAD*0600'));
        % Know a priori what the bandwidth of the coefficients is
        Ldata=180;
    end
elseif  strcmp(Pcenters,'CSR')
    if strcmp(Rlevel,'RL04')
        datanames=ls2cell(fullfile(ddir1,'GAD*0060_0004'));
    elseif strcmp(Rlevel,'RL05')
        datanames=ls2cell(fullfile(ddir1,'GSM*0005'));
        % Know a priori what the bandwidth of the coefficients is
        Ldata=100;
    elseif strcmp(Rlevel,'RL06')
        % 5/12/2022 GR-FO only now version 6.1
        datanames=ls2cell(fullfile(ddir1,'GAD*0600'));
        % Know a priori what the bandwidth of the coefficients is
        Ldata=180;
    end
elseif  strcmp(Pcenters,'JPL')
    if strcmp(Rlevel,'RL05')
        datanames=ls2cell(fullfile(ddir1,'GAD*0005'));
        % Know a priori what the bandwidth of the coefficients is
        Ldata=100;
    elseif strcmp(Rlevel,'RL06')
        % 5/12/2022 GR-FO only now version 6.1
        datanames=ls2cell(fullfile(ddir1,'GAD*0600'));
        % Know a priori what the bandwidth of the coefficients is
        Ldata=180;
    end
end

% Initialize
nmonths = length(datanames);
thedates = zeros(1,nmonths);

fid=fopen(ouput_filename,'wt');

fprintf(fid,'An Additional File referencing to GRACE Technical Note 13\n');
fprintf(fid,'\n');
fprintf(fid,'TITLE:\n');
fprintf(fid,'\tMonthly degree-1 (geocenter) gravity coefficients from GAD product, generated from\n');
fprintf(fid,'\tGRACE (01-2003 - 12/2024) %s solutions.\n',Rlevel);
fprintf(fid,'\tLast reported data point: %s.\n',datetime('today'));
fprintf(fid,'\n');
fprintf(fid,'DESCRIPTION:\n');
fprintf(fid,'\tThis file contains estimates of degree-1 gravity coefficients for the GRACE Release-%s data\n',Rlevel(3:4));
fprintf(fid,'\tbased on GAD products for ocean application.\n');
fprintf(fid,'\n');
fprintf(fid,'\tReference in GRACE Technical Note 13:\n');
fprintf(fid,'\t\tFor ocean applications, add the degree 1 values found in the ocean bottom pressure product \n');
fprintf(fid,'\t\t(GAD) to these values. If you are interested in the full ocean, land, atmosphere geocenter, \n');
fprintf(fid,'\t\tadd the degree 1 values of the atmosphere/ocean product (GAC) to these values. \n');
fprintf(fid,'\n');
fprintf(fid,'SPECIAL NOTES:\n');
fprintf(fid,'\tThe anomalous signals are relative to the mean field from 2003-2014.\n');
fprintf(fid,'\n');
fprintf(fid,'FORMAT:\n');
fprintf(fid,'\tThe data format is similar to that for GSM, GAD, and GAC files distributed by PODAAC, see\n');
fprintf(fid,'\tGRACE Level-2 handbook for more information.\n');
fprintf(fid,'\n');
fprintf(fid,'# EGM Coefficient Record 2 \n');
fprintf(fid,'  variables:\n');
fprintf(fid,'    record_key:\n');
fprintf(fid,'      key_name              : GRCOF2 \n');
fprintf(fid,'      long_name             : Earth Gravity Spherical Harmonic Model Format Type\n');
fprintf(fid,'      data_type             : string\n');
fprintf(fid,'      comment               : 1st column\n');
fprintf(fid,'    degree_index:\n');
fprintf(fid,'      key_name              : spherical harmonic degree l \n');
fprintf(fid,'      data_type             : integer\n');
fprintf(fid,'      comment               : 2nd column\n');
fprintf(fid,'    order_index:\n');
fprintf(fid,'      key_name              : spherical harmonic degree l \n');
fprintf(fid,'      data_type             : integer\n');
fprintf(fid,'      comment               : 3rd column\n');
fprintf(fid,'    clm:\n');
fprintf(fid,'      key_name              : Clm coefficient; cosine coefficient for degree l and order m\n');
fprintf(fid,'      data_type             : double precision\n');
fprintf(fid,'      comment               : 4th column\n');
fprintf(fid,'    slm:\n');
fprintf(fid,'      key_name              : Slm coefficient; cosine coefficient for degree l and order m\n');
fprintf(fid,'      data_type             : double precision\n');
fprintf(fid,'      comment               : 5th column\n');
fprintf(fid,'    clm_std_dev: \n');
fprintf(fid,'      key_name              : standard deviation of Clm\n');
fprintf(fid,'      data_type             : double precision\n');
fprintf(fid,'      comment               : 6th column\n');
fprintf(fid,'    slm_std_dev: \n');
fprintf(fid,'      key_name              : standard deviation of Slm\n');
fprintf(fid,'      data_type             : double precision\n');
fprintf(fid,'      comment               : 7th column\n');
fprintf(fid,'    epoch_begin_time: \n');
fprintf(fid,'      long_name             : epoch begin of Clm, Slm coefficients\n');
fprintf(fid,'      time_format           : yyyymmdd.hhmm\n');
fprintf(fid,'      data_type             : double precision\n');
fprintf(fid,'      comment               : 8th column\n');
fprintf(fid,'    epoch_stop_time: \n');
fprintf(fid,'      long_name             : epoch stop of Clm, Slm coefficients\n');
fprintf(fid,'      time_format           : yyyymmdd.hhmm\n');
fprintf(fid,'      data_type             : double precision\n');
fprintf(fid,'      comment               : 9th column\n');
fprintf(fid,'\n');
fprintf(fid,'end of header ===============================================================================\n');

% Loop over the months
for index = 1:nmonths
    disp(['Processing month number: ' num2str(index)]);

    % load geopotential coefficients
    fname1=fullfile(ddir1,datanames{index});

    % Open and scan the file (data from all three centers is 10 columns)
    fid2 = fopen(fname1);
    C = textscan(fid2,'%s%s%s%s%s%s%s%s%s%s');
    fclose(fid2);

    % Only grab the lines for GRCOF2
    Carray = cat(3,C{:});
    I = strmatch('GRCOF2',Carray(:,1,1),'exact');
    Carray = squeeze(Carray(I,1,:));

    % columns 2-9, and as format double
    Carray = Carray(:,2:9);
    lmcosi_month=cellfun(@str2num,Carray);
    % This should be addmup(Ldata)
    [m,n] = size(lmcosi_month);

    d1o0=find(lmcosi_month(:,1)==1 & lmcosi_month(:,2)==0);
    d1o1=find(lmcosi_month(:,1)==1 & lmcosi_month(:,2)==1);
    fprintf(fid,'GRCOF2\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',Carray{d1o0,1},...
        Carray{d1o0,2},Carray{d1o0,3},Carray{d1o0,4},Carray{d1o0,5},...
        Carray{d1o0,6},Carray{d1o0,7},Carray{d1o0,8});
    fprintf(fid,'GRCOF2\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',Carray{d1o1,1},...
        Carray{d1o1,2},Carray{d1o1,3},Carray{d1o1,4},Carray{d1o1,5},...
        Carray{d1o1,6},Carray{d1o1,7},Carray{d1o1,8});
end
fclose(fid);

;
end
