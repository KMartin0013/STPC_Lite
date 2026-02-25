function maybe_plot_region(XY_buffer, XY_ori, FIG_Attach, ...
    ll_buffer, figDir, doPlot)
% MAYBE_PLOT_REGION  Optionally plot buffered and original basin contours.

    if ~doPlot
        return;
    end

    f0 = figure('Visible','on');
    geoshow(XY_buffer(:,2), XY_buffer(:,1), 'DisplayType', 'polygon', ...
        'FaceColor', 'yellow'); hold on;

    if XY_ori(1,1) > 180
        plot(XY_ori(:,1)-360, XY_ori(:,2), '--b', 'LineWidth', 3);
    else
        plot(XY_ori(:,1), XY_ori(:,2), '--b', 'LineWidth', 3);
    end
    title(['Study area and buffer zone (' num2str(ll_buffer) char(176) ')']);

    tif_name0 = ['Fig0_' FIG_Attach '.tif'];
    if ~exist(figDir, 'dir'); mkdir(figDir); end
    print(f0,'-dtiff','-r125', fullfile(figDir, tif_name0));
%     close(f0);
end