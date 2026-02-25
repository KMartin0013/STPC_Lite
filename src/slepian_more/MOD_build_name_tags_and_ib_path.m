function [Dataproduct, FIG_Attach_Mod] = MOD_build_name_tags_and_ib_path( ...
    ii, oo, rr, Area, Lwindow);
% BUILD_NAME_TAGS_AND_IB_PATH  Build figure/file name prefix and IB file path.

    if isempty(rr)
        Dataproduct = {num2str(ii),num2str(oo)};
        FIG_Attach_Mod=sprintf('I%sO%s_%s_%s',...
            Dataproduct{1},Dataproduct{2},Area,num2str(Lwindow));
        
    elseif length(rr)==1
        Dataproduct = {num2str(ii),num2str(oo),num2str(rr)};
        FIG_Attach_Mod=sprintf('I%sO%sR%s_%s_%s',...
            Dataproduct{1},Dataproduct{2},Dataproduct{3},Area,num2str(Lwindow));
    end

end