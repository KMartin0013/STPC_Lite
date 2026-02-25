function setup_paths(ifilesRoot)
% SETUP_PATHS  Add dependent code directories to MATLAB path
%
% Usage:
%   setup_paths(ifilesRoot)
%
% Input:
%   ifilesRoot - char, path to the original IFILES root directory

    arguments
        ifilesRoot (1,:) char
    end


    create_paths(fullfile(ifilesRoot,'COASTS'));
    create_paths(fullfile(ifilesRoot,'GLMALPHA'));
    create_paths(fullfile(ifilesRoot,'KERNELC'));
    create_paths(fullfile(ifilesRoot,'LEGENDRE'));

    try 
    % https://github.com/strawpants/GRACE-filter (DDK software)
    check_softwares(fullfile(ifilesRoot,'src','required_softwares','GRACE-filter-master'));
    catch

    warning('You do not download the DDK procedure from ''https://github.com/strawpants/GRACE-filter''. Some smoothing procedure may not work.')
    end

    % https://geoweb.princeton.edu/people/simons/software.html (Slepian softwares)
    check_softwares(fullfile(ifilesRoot,'src','required_softwares','slepian_alpha-master'));
    check_softwares(fullfile(ifilesRoot,'src','required_softwares','slepian_bravo-master'));
    check_softwares(fullfile(ifilesRoot,'src','required_softwares','slepian_delta-master'));
    check_softwares(fullfile(ifilesRoot,'src','required_softwares','slepian_zero-master'));

end