function plot_mssa_evalues_and_turning_points(mssa_results, V, Max_S, S_bou, ...
    N_bou, turn_V_four, turn_MSSA_four, MSSA_evalues_sumup, ...
    Turning_number, used_sta_SN, Attach);

intit_num           = mssa_results.intit_num;
S                   = mssa_results.S;
M_rec               = mssa_results.M_rec;

MaxNumChanges_MSSA=Turning_number; % how many turning point do you want?
MaxNumChanges_V=Turning_number; % how many turning point do you want?

option_sta=["mean","rms","std","linear"];

used_sta_V=used_sta_SN; % what kind of turning point do you want? (default is RSS)sed_sta_V=4; % what kind of turning point do you want? (default is RSS)
used_sta_MSSA=used_sta_SN; % what kind of turning point do you want? (default is RSS)

%% combination of V and MSSA
font_Size=14;F_Position=[100,250,1400,500];
F_gap=[.03 .08];F_marg_h=[.12 .05];F_marg_w=[.07 .02];

f1c=figure

set(gcf,'Position',F_Position)

ha=tight_subplot(1,2,F_gap,F_marg_h,F_marg_w);

axes(ha(1))

if S_bou~=0
    disp(['You choose the boundary restriction for S: ' num2str(S_bou)])
    he_gap1=area([0,Max_S],[S_bou 1-2*S_bou S_bou;S_bou 1-2*S_bou S_bou],0);
    he_gap1(2).FaceColor = [1,1,1];
    he_gap1(2).EdgeColor = [1,1,1];

    he_gap1(1).FaceColor = [0.4,0.4,0.4];
    he_gap1(1).EdgeColor = [0.4,0.4,0.4];
    he_gap1(1).FaceAlpha = 0.4;
    he_gap1(1).EdgeAlpha = 0.7;
    he_gap1(3).FaceColor = [0.4,0.4,0.4];
    he_gap1(3).EdgeColor = [0.4,0.4,0.4];
    he_gap1(3).FaceAlpha = 0.4;
    he_gap1(3).EdgeAlpha = 0.7;

    hold on
end

p1=plot(1:Max_S,V(1:Max_S),'Color',[0 0.4470 0.7410],'LineWidth',2);
hold on
for i=1:MaxNumChanges_V
    text(turn_V_four(i,used_sta_V)+2,V(turn_V_four(i,used_sta_V)),['S(' num2str(i) ')'],...
        'FontSize',font_Size,'Color','r')
    hold on
end

p2=plot([0 turn_V_four(1,used_sta_V)],[mean(V(1:turn_V_four(1,used_sta_V))) mean(V(1:turn_V_four(1,used_sta_V)))],...
    'Color','k','LineWidth',1.5,'LineStyle','--');
hold on
for i=1:MaxNumChanges_V-1
    plot([turn_V_four(i,used_sta_V) turn_V_four(i,used_sta_V)],[0 1],...
        'Color',[0.8 0.8 0.8],'LineWidth',1)
    hold on
    plot([turn_V_four(i,used_sta_V) turn_V_four(i+1,used_sta_V)],[mean(V(turn_V_four(i,used_sta_V):turn_V_four(i+1,used_sta_V))) mean(V(turn_V_four(i,used_sta_V):turn_V_four(i+1,used_sta_V)))],...
        'Color','k','LineWidth',1.5,'LineStyle','--');
    hold on
end
plot([turn_V_four(MaxNumChanges_V,used_sta_V) turn_V_four(MaxNumChanges_V,used_sta_V)],[0 1],...
    'Color',[0.8 0.8 0.8],'LineWidth',1)
hold on
plot([turn_V_four(MaxNumChanges_V,used_sta_V) length(V)],[mean(V(turn_V_four(MaxNumChanges_V,used_sta_V):length(V))) mean(V(turn_V_four(MaxNumChanges_V,used_sta_V):length(V)))],...
    'Color','k','LineWidth',1.5,'LineStyle','--');
hold on

p3=plot(turn_V_four(:,used_sta_V),V(turn_V_four(:,used_sta_V)),'o','MarkerSize',5,...
    'MarkerFaceColor','r','MarkerEdgeColor','k');
hold on

% xlim([0,Max_S])

ylabel({'\bf Concentration ratio ($\lambda\left(\alpha\right)$)','\bf of Slepian conversion'},...
    'Interpreter','latex','FontWeight','bold')
xlabel('$\alpha$','Interpreter','latex','FontWeight','bold')
legend([p3],{['S(1~' num2str(MaxNumChanges_V) ')']}, ...
    'NumColumns',2,'Position',[0.3,0.7,0.2,0.1],'FontSize',font_Size+2,'box','off')

set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

% set(gca, 'xticklabel', get(gca, 'xtick'))
xlim([0,Max_S])

axes(ha(2))
% f1b=figure

max_mode=min(150,intit_num*M_rec);

if N_bou~=0
    disp(['You choose the boundary restriction for N: ' num2str(N_bou)])
    he_gap1=area([0,max_mode],[N_bou 1-2*N_bou N_bou;N_bou 1-2*N_bou N_bou],0);
    he_gap1(2).FaceColor = [1,1,1];
    he_gap1(2).EdgeColor = [1,1,1];

    he_gap1(1).FaceColor = [0.4,0.4,0.4];
    he_gap1(1).EdgeColor = [0.4,0.4,0.4];
    he_gap1(1).FaceAlpha = 0.4;
    he_gap1(1).EdgeAlpha = 0.7;
    he_gap1(3).FaceColor = [0.4,0.4,0.4];
    he_gap1(3).EdgeColor = [0.4,0.4,0.4];
    he_gap1(3).FaceAlpha = 0.4;
    he_gap1(3).EdgeAlpha = 0.7;

    hold on
end

colorr=colormap('winter');
dd=linspace(1,256,S);
map_colori=colorr(floor(flipud(dd')),:);TickLabels_c={};
for j=1:S
    %     color_line=repmat(0.2+0.02*j,1,3);
    color_line=map_colori(j,:);
    plot(1:max_mode,MSSA_evalues_sumup(1:max_mode,j),'color',color_line,'LineWidth',1);
    hold on

    p3=plot(turn_MSSA_four{j}(:,used_sta_MSSA),MSSA_evalues_sumup(turn_MSSA_four{j}(:,used_sta_MSSA),j),'o','MarkerSize',5,...
        'MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerEdgeColor','k');
    hold on

    xlim([0,max_mode])

    for i=1:MaxNumChanges_MSSA
        turning_MSSA(j,i)=turn_MSSA_four{j}(i,used_sta_MSSA);
        turning_MSSA_evalues(j,i)=MSSA_evalues_sumup(turn_MSSA_four{j}(i,used_sta_MSSA),j);
    end

    TickLabels_c{j}=num2str(j);
end

colormap(gca,map_colori)
%         colorbar('Ticks',linspace(0,1,iter_num(i)),...
%             'TickLabels',TickLabels_c)
cb=colorbar('Ticks',[1/S/2:10/S:1,1],...
    'TickLabels',TickLabels_c([1:10:end,end]))

for i=1:MaxNumChanges_MSSA
    plot(turning_MSSA(1:S,:),turning_MSSA_evalues(1:S,:),'k')

    hold on
    text(max(turning_MSSA(:,i))+1,min(turning_MSSA_evalues(:,i))-0.01,['N(' num2str(i) ')'],...
        'FontSize',font_Size,'Color',[0.9290 0.6940 0.1250])
    hold on
end

% cb = colorbar;
cb.Position = [0.9 0.2 0.02,0.3];
title(cb,'$\alpha$','Interpreter','latex','fontsize',14);

% xlim([0,max_mode])

ylabel({'\bf Sum of temporal eigenvalues ($\sum_{\beta}{\widetilde{\lambda}}_\alpha(\beta)$)', ...
    '\bf of MSSA'},'Interpreter','latex','FontWeight','bold')
xlabel('$\beta$','Interpreter','latex','FontWeight','bold')

legend([p3],{['N(1~' num2str(MaxNumChanges_MSSA) ')']}, ...
    'NumColumns',2,'Position',[0.8,0.7,0.2,0.1],'FontSize',font_Size+2,'box','off')

set(gca,'FontName','Times New Roman','FontSize',font_Size,'TickLength',[0.004,0.035]);

xlim([0,max_mode])
ylim([0,1])
% set(gca, 'xticklabel', get(gca, 'xtick'))
% xlim([0,Max_SSF])

% set(gcf,'Position',F_Position)

tif_name1c=['TurningPoint_V_MSSA.tif'];
print(f1c,'-dtiff','-r300',fullfile(Attach.fig_path_ALL,tif_name1c));

% save(fullfile(getenv('IFILES'),['Results_' parrent_path],['Main_' Attach_ALL ...
%     '_TP' num2str(Turning_number) '_M' num2str(M_rec) '_Breakingpoint.mat']),...
%     'MaxNumChanges_MSSA','MaxNumChanges_V','turn_V_four','turn_MSSA_four',...
%     'N_gap','M_gap','M_rec','Turning_number','RC_RMS_Seq','V',...
%     'used_sta_V','used_sta_MSSA','MSSA_evalues_sumup','Max_S','test_M_rec',...
%     'threshold_S','threshold_N')

% figure
% plot(XY(:,1),XY(:,2))
% hold on
% plot(XY_ori(:,1),XY_ori(:,2))

disp('Now we have already determined the breaking points of S and N. ')
disp('Then, we need the remove the noise by statistical method at different significance level p.')

end
