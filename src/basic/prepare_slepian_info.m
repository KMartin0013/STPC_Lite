function slepianInfo = prepare_slepian_info(config)
% PREPARE_SLEPIAN_INFO  Build struct to be saved as Slepian_Information.mat
%
% Usage:
%   slepianInfo = prepare_slepian_info(config)

    slepianInfo = struct();
    
    
        % Do you have already determine the buffer zone?
    if isnumeric(config.buffer_deg)

        % Yes. Just use it.
        slepianInfo.buffer_deg=config.buffer_deg;

    elseif strcmp(config.buffer_deg,'Auto')
        % NO. We need to determine it. 

        slepianInfo.buffer_deg='Auto';

        % Do you have a preferred buffer zone range?
        if isfield(config,'groupBuffer') 

            if numel(config.groupBuffer)==1

                error('You should provide more than one preffered buffer zones.')
            elseif all(config.groupBuffer >= 0) || all(config.groupBuffer <= 0) 

                warning('We will find the suitable buffer zone based on your preffered ranges.')
            
                slepianInfo.group_buffer = config.groupBuffer;
            else

                error(['You should provide consistent signs (positive' ...
                    ' or negative) of preffered buffer zones.'])
            end

        else
            warning('We will find the suitable buffer zone based on empirical ranges.')

            if strcmp(config.landOrOcean,'land') || strcmp(config.landOrOcean,'ice')

                slepianInfo.group_buffer = 0:0.5:1.5;
            elseif strcmp(config.landOrOcean,'ocean')

                slepianInfo.group_buffer = 0:-0.5:-1.5;
            end
        end

    else

        error('You should provide a buffer zone OR set it as ''Auto'' and provide preffered buffer zones ranges in ''config.groupBuffer''. ')
    end

    
    % Other parameters
    slepianInfo.buffer_deg         = config.buffer_deg;
    slepianInfo.S_choice           = config.S_choice;
    slepianInfo.Radius             = config.Radius;
    slepianInfo.phi                = config.phi;
    slepianInfo.theta              = config.theta;
    slepianInfo.omega              = config.omega;
    slepianInfo.artificial_months  = config.artificialMonths;
end