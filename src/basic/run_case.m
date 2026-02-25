function STPC_results = run_case(config)
% RUN_CASE_YANGTZE  Main controller for the cases
%
% Usage:
%   results = run_case(config)
%
% Input:
%   config - struct containing paths, region, Slepian, MSSA, etc.
%
% Output:
%   results - struct collecting results of each stage
%             (or you can rely on files written to disk only)

    arguments
        config struct
    end

    %---------------------------%
    % 1. Paths and output dirs  %
    %---------------------------%
    setup_paths(config.ifilesRoot);

    % Make sure output directories exist
    if ~exist(config.resultDir, 'dir')
        mkdir(config.resultDir);
    end
    if ~exist(config.figureDir, 'dir')
        mkdir(config.figureDir);
    end

    % Derived field: full codeVersion including basin tag
    codeVersion = [config.codeVersionBase config.TH_ori];

    %---------------------------%
    % 2. Save Basic information %
    %---------------------------%
    basicInfo = prepare_basic_info(codeVersion, config);
    save(fullfile(config.resultDir,'Basic_Information.mat'), ...
        '-struct', 'basicInfo');

    %---------------------------------------------%
    % 3. Prepare and run Slepian info procedure   %
    %---------------------------------------------%
    
    num_givIns = numel(config.givInstitu);

    slepianInfo = prepare_slepian_info(config);

    if strcmp(slepianInfo.buffer_deg,'Auto')

        fprintf(['==== Running Slepian procedure (det_opt_buffer) to determine buffer zone. ====\n']);

        buffer_deg=det_opt_buffer(basicInfo, slepianInfo);

        slepianInfo.buffer_deg=buffer_deg;

        fprintf('==== End (buffer zone: %s) ====\n', num2str(slepianInfo.buffer_deg));

    end

    % buffer string used in file/figure names
    if slepianInfo.buffer_deg >= 0
        buffer_str = num2str(slepianInfo.buffer_deg);
        buffer_str(buffer_str=='.') = 'p';
    else
        buffer_str = ['neg' num2str(abs(slepianInfo.buffer_deg))];
        buffer_str(buffer_str=='.') = 'p';
    end

    slepianInfo.buffer_str = buffer_str;

    save(fullfile(config.resultDir,'Slepian_Information.mat'),...
        '-struct', 'slepianInfo');

%     slepian_results=struct();
    for Ins = 1 : num_givIns

        dataProduct   = {char(config.givInstitu(Ins)),config.Institu_ver, ...
            config.Lwindow,config.landOrOcean,config.GIA};

        fprintf(['==== Running Slepian procedure (run_slepian_main) for ''%s'' ====\n'] , ...
            char(config.givInstitu(Ins)));
        slepian_result=run_slepian_main(basicInfo, slepianInfo, dataProduct);

        disp(['SLEPIAN: The time series ''Code'' and ''Area-weighted'' should be ', ...
            'generally consistent. Otherwise, check c11cmn & buffer zones.']);
%         pause;  % remove this line if you do not want an interactive pause

        fprintf(['==== End for ''%s''' ' ====\n'], char(config.givInstitu(Ins)))

        slepian_results(Ins) = slepian_result;
    end

    %---------------------------%
    % 5. MSSA settings & run    %
    %---------------------------%

    mssaInfo = prepare_mssa_info(config);
    save(fullfile(config.resultDir,'MSSA_Information.mat'), ...
        '-struct', 'mssaInfo');

    mssa_smooth_results = run_mssa_smooth(basicInfo, slepianInfo, mssaInfo, ...
        slepian_results);  % internally wraps Lite_MSSA_V1

    mssa_results = run_mssa_main(basicInfo, slepianInfo, mssaInfo, ...
        slepian_results);  % internally wraps Lite_MSSA_V1

    %---------------------------%
    % 6. STPC filtering         %
    %---------------------------%
    STPCInfo = prepare_stpc_info(config);
    save(fullfile(config.resultDir,'STPC_Information.mat'),...
        '-struct', 'STPCInfo');

    STPC_p_results = run_stpc_main(basicInfo, STPCInfo, slepian_results, ...
        mssa_results);  % internally wraps Lite_MSSA_V1

    STPC_results = [STPC_p_results, mssa_smooth_results];

end