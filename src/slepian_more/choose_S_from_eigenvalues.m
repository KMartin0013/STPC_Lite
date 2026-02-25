function [S_shannon, S, S_sig] = ...
    choose_S_from_eigenvalues(S_choice, Max_S, V, Lwindow, XY_buffer)
% CHOOSE_S_FROM_EIGENVALUES  Determine S and S_sig from eigenvalue spectrum.
%
% Inputs:
%   S_choice   - integer flag specifying strategy
%   Max_S      - maximum S when S_choice == 6
%   V          - eigenvalues of Slepian basis
%   Lwindow    - maximum spherical harmonic degree
%   XY_buffer  - buffered region polygon
%   FIG_Attach - base name tag to be extended
%
% Outputs:
%   S_shannon  - Shannon number
%   S          - chosen number of Slepian functions
%   S_sig      - either S or [N1, N2] (for smoothing cases)
%   FIG_Attach - updated tag including S information

    S_shannon = round((Lwindow + 1)^2 * spharea(XY_buffer));
    S1 = []; S2 = [];

    switch S_choice
        case 0
            % Use Shannon number directly
            S     = S_shannon;
            S_sig = S;
        case 1
            % Eigenvalues > 0.1
            S     = sum(V > 0.1);
            S_sig = S;
        case 2
            % Provide [N1, N2] where N1 = Shannon, N2 = V>0.1
            S  = sum(V > 0.1);
            S1 = S_shannon;
            S2 = S;
            S_sig = [S1, S2];
        case 3
            % Smooth all Slepian functions
            S  = S_shannon;
            S1 = 0;
            S2 = S;
            S_sig = [S1, S2];
        case 4
            % Eigenvalues > 0.3
            S     = sum(V > 0.3);
            S_sig = S;
        case 5
            % Eigenvalues > 0.01
            S     = sum(V > 0.01);
            S_sig = S;
        case 6
            % Use Max_S as S
            S     = Max_S;
            S_sig = S;
        otherwise
            warning('Unknown S_choice. Using Shannon number.');
            S     = S_shannon;
            S_sig = S;
    end

    % Extend name tag with S info
%     if numel(S_sig) > 1
%         % If you provide two elements in 'S_sig', there should be additional
%         % Gaussian smoothing on high-order SSF, see Ref2 for details
%         FIG_Attach = [FIG_Attach '_S' num2str(S_sig(1)) 't' num2str(S_sig(2))];
%     else
%         
%         FIG_Attach = [FIG_Attach '_S' num2str(S)];
%     end

end