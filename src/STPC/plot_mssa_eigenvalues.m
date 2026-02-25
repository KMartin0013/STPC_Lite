function plot_mssa_eigenvalues( ...
    mssa_Sort_coff_sn,intit_num,S_rec,V, ...
    TP_ss,TP_nn,Turning_number,M_rec,Noise_SigLev,Attach)
% -------------------------------------------------------------------------
% plot_MSSA_eigenvalues
%
% This function reproduces the two plotting blocks:
%   (1) "MSSA_Eigenvalues_box_Frequency_V_TP..."  (figure f1a)
%   (2) "MSSA_Eigenvalues_box_noise_Frequency_V_TP..." (figure f1b)
%
% INPUT:
%   mssa_Sort_coff_sn : struct array containing MSSA coefficient evaluation
%                       results. Must contain field:
%                       - fillM_reconcoffs_eva (eigenvalue-related matrix)
%   intit_num         : MSSA iteration index (will use intit_num+1 entry)
%   S_rec             : number of Slepian coefficients used (S in script)
%   V                 : vector of Slepian concentration ratios (lambda)
%                       length(V) >= S_rec
%   TP_ss             : truncation number S used in title
%   TP_nn             : cutoff number N used in title
%   Turning_number    : turning number used in output file naming
%   M_rec             : embedding dimension / reconstruction parameter
%                       used in output file naming
%   Noise_SigLev      : noise significance level (0–1), used in text box
%                       and output file naming
%   Attach            : directory where figures (tif) will be saved
%
% OUTPUT:
%   f1a, f1b          : figure handles to the two created figures
%
% NOTE:
%   - This function uses external utilities: tight_subplot, hatchfill2.
%   - All plotting parameters and variable names follow the original code.
%   - No computation logic was changed, only encapsulated as a function.
% -------------------------------------------------------------------------

%% Common plotting parameters

font_Size = 16;
MarkSize = 2; %#ok<NASGU>
F_Position = [100,50,700,500];
abc = 'abcdefghijklmnopqrstuvwxyz'; %#ok<NASGU>

% Blind-friendly color palette
colorr_blind = [213,94,0; ...
                204,121,167; ...
                0,158,115; ...
                230,159,0; ...
                213,94,0; ...
                240,228,66; ...
                0,114,178; ...
                86,180,233]/255;

ybou1 = [0,1.19];
F_gap = [.02 .03];
F_marg_h = [.12 .08];
F_marg_w = [.10 .10];

% eigenvalue-related evaluation for plotting
used_eva = mssa_Sort_coff_sn(intit_num+1).fillM_reconcoffs_eva;

%% ------------------------------------------------------------------------
%  Figure f1a: eigenvalues distribution + Slepian concentration ratio
% -------------------------------------------------------------------------
f1a = figure;
left_color  = [0 0 0];
right_color = [255 0 0]/255;
set(f1a,'defaultAxesColorOrder',[left_color;right_color])

ha = tight_subplot(1,1,F_gap,F_marg_h,F_marg_w);
axes(ha(1));

yyaxis left
x = []; %#ok<NASGU>

% We do not want IB to show here; S_ol-related code is omitted
X = 1:S_rec;
% Y columns: [long-term, annual, short-term, noise]
Y = [used_eva(1:S_rec,[3,4,6]), used_eva(1:S_rec,1)-used_eva(1:S_rec,2)];

bb = bar(X,Y,'stacked');
set(bb(1),'FaceColor',colorr_blind(1,:));  % long-term
set(bb(2),'FaceColor',colorr_blind(3,:));  % annual
set(bb(3),'FaceColor',colorr_blind(4,:));  % short-term
set(bb(4),'FaceColor',[169,169,169]/255);  % noise

ylim(ybou1);
xlim([0,S_rec+1])

ylabel('sum of normalized eigenvalues in M-SSA ($\sum\widetilde{\lambda}$)', ...
       'Interpreter','latex')
xlabel('Spherical slepian coefficients')

yyaxis right
y1 = plot(1:S_rec,V(1:S_rec),'r','LineWidth',2); %#ok<NASGU>

hold on
% Shannon-related annotation was commented in original code

ylim(ybou1);
xlim([0,S_rec+1])
grid on

ylabel('concentration ratio of spherical Slepian functions (\lambda)')

legend({'long-term','annual','short-term','noise','\lambda'}, ...
    'NumColumns',3,'Position',[0.4,0.84,0.4,0.03],'box','off')

title(['TP-S(' num2str(TP_ss) ') and TP-N(' num2str(TP_nn) ')'])

set(gca,'FontName','Times New Roman', ...
        'FontSize',font_Size-2, ...
        'TickLength',[0.004,0.035]);

xticks(get(gca,'xtick'))
set(gca,'xticklabel',get(gca,'xtick'))

set(gcf,'Position',F_Position)

tif_name1a = ['MSSA_Eigenvalues_box_Frequency_V_TP' num2str(Turning_number) ...
    '_S' num2str(TP_ss) '_N' num2str(TP_nn) '_M' num2str(M_rec) ...
    '_Sig' num2str(Noise_SigLev) '.tif'];

print(f1a,'-dtiff','-r300',fullfile(Attach.fig_path_ALL,tif_name1a));

pause(0.5)
close(f1a)
%% ------------------------------------------------------------------------
%  Figure f1b: eigenvalues distribution including noise frequencies
% -------------------------------------------------------------------------

% Re-use the same common settings
used_eva = mssa_Sort_coff_sn(intit_num+1).fillM_reconcoffs_eva;

f1b = figure;
left_color  = [0 0 0];
right_color = [255 0 0]/255;
set(f1b,'defaultAxesColorOrder',[left_color;right_color])

ha = tight_subplot(1,1,F_gap,F_marg_h,F_marg_w);
axes(ha(1));

x = []; %#ok<NASGU>

X = 1:S_rec;
% Y columns:
%   [3,4,6]  : long-term, annual, short-term
%   [8,9,11]: long-term-noise, annual-noise, short-term-noise
%   last    : residual part (1 - total explained)
Y = [used_eva(1:S_rec,[3,4,6]), ...
     used_eva(1:S_rec,[8,9,11]), ...
     ones(size(used_eva(1:S_rec,1))) - used_eva(1:S_rec,1)];

bb = bar(X,Y,'stacked');
set(bb(1),'FaceColor',colorr_blind(1,:));  % long-term
set(bb(2),'FaceColor',colorr_blind(3,:));  % annual
set(bb(3),'FaceColor',colorr_blind(4,:));  % short-term

% hatch patterns for noise-related components (4:6)
hatchfill2(bb(4),'single','HatchAngle',60, ...
    'hatchcolor',colorr_blind(1,:),'HatchLineWidth',1.5);
hatchfill2(bb(5),'single','HatchAngle',60, ...
    'hatchcolor',colorr_blind(3,:),'HatchLineWidth',1.5);
hatchfill2(bb(6),'single','HatchAngle',60, ...
    'hatchcolor',colorr_blind(4,:),'HatchLineWidth',1.5);

for b = 4:6
    bb(b).FaceColor = 'none';
end

% residual (noise) as solid gray
set(bb(7),'FaceColor',[169,169,169]/255);

ylim(ybou1);
xlim([0,S_rec+1])

grid on

ylabel('sum of normalized eigenvalues in M-SSA ($\sum\widetilde{\lambda}$)', ...
       'Interpreter','latex')
xlabel('Spherical slepian coefficients')

% Right axis with lambda was commented in original code

legend([bb(1),bb(2),bb(3),bb(7)], ...
       {'long-term','annual','short-term','residual'}, ...
       'NumColumns',3,'Location','northeast','box','off')

tt = text(1,1.1,['\gamma = ' num2str((1-Noise_SigLev)*100) '%'], ...
    'FontSize',font_Size-2, ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','left', ...
    'BackgroundColor','white', ...
    'Color','r', ...
    'EdgeColor','black'); %#ok<NASGU>

title({['Truncation number S(' num2str(TP_ss) ...
        ') and cutoff number N(' num2str(TP_nn) ')']}, ...
      'FontSize',font_Size-3)

set(gca,'FontName','Times New Roman', ...
        'FontSize',font_Size-2, ...
        'TickLength',[0.004,0.035]);

xticks(get(gca,'xtick'))
set(gca,'xticklabel',get(gca,'xtick'))

set(gcf,'Position',F_Position)

tif_name1b = ['MSSA_Eigenvalues_box_noise_Frequency_V_TP' num2str(Turning_number) ...
    '_S' num2str(TP_ss) '_N' num2str(TP_nn) '_M' num2str(M_rec) ...
    '_Sig' num2str(Noise_SigLev) '.tif'];

print(f1b,'-dtiff','-r300',fullfile(Attach.fig_path_ALL,tif_name1b));

pause(0.5)
close(f1b)
end
