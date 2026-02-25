function basicInfo = prepare_basic_info(codeVersion, config)
% PREPARE_BASIC_INFO  Build struct to be saved as Basic_Information.mat
%
% Usage:
%   basicInfo = prepare_basic_info(codeVersion, config)

    basicInfo = struct();
    basicInfo.Code_Version  = codeVersion;
    basicInfo.TH_ori        = config.TH_ori;
    basicInfo.Area          = config.areaName;
    basicInfo.Lwindow       = config.Lwindow;
    basicInfo.land_or_ocean = config.landOrOcean;
    basicInfo.c11cmn        = config.c11cmn;
    basicInfo.Max_S         = config.Max_S;
    basicInfo.ddir1         = config.resultDir;
    basicInfo.ddir2         = config.figureDir;
    basicInfo.ifilesRoot    = config.ifilesRoot;
    basicInfo.redo          = config.redo;
    basicInfo.plotProcess   = config.plotProcess;
    basicInfo.saveAddData   = config.saveAddData;
    basicInfo.note          = config.note;
    basicInfo.Smooth        = config.Smooth;

    basicInfo.Lwindow       = config.Lwindow;
    basicInfo.GIA           = config.GIA;
    basicInfo.Institu_ver   = config.Institu_ver;
    basicInfo.use_institu   = config.use_institu;

end