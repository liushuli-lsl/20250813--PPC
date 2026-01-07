%% compare_x1_to_x5.m
clc;
close all;

num_cases = 1;
% tp_list = [1, 3, 5];
% num_cases = numel(tp_list);
colors = lines(7);
line_styles = {'-', '--', '-.', ':','-', '--','-.', ':',}; % 不同线型

% 准备存储
all_tspan = cell(num_cases,1);
all_e_q1   = cell(num_cases,1);
all_e_q2   = cell(num_cases,1);
all_e_dq1   = cell(num_cases,1);
all_e_dq2   = cell(num_cases,1);
all_tau1   = cell(num_cases,1);
all_tau2   = cell(num_cases,1);
all_u1  = cell(num_cases,1);
    all_u2  = cell(num_cases,1);
all_q1     = cell(num_cases,1);
all_q2     = cell(num_cases,1);
all_qd1    = cell(num_cases,1);
all_qd2    = cell(num_cases,1);
all_dq1     = cell(num_cases,1);
all_dq2     = cell(num_cases,1);
all_dqd1    = cell(num_cases,1);
all_dqd2    = cell(num_cases,1);
all_rho1    = cell(num_cases,1);
all_rho2    = cell(num_cases,1);
for k = 1:num_cases
     fname = sprintf('x%d.mat', 3);
    S = load(fname, ...
        'tspan','e_q','e_dq','tau_mat','u_mat','q_use','qd_mat','dq_use','dqd_mat','rho_mat');
    t_full = S.tspan(1:end-1);
idx    = t_full <= 20;        % 逻辑索引
    all_tspan{k} = S.tspan(idx);           % 1×N
    all_e_q1{k}  = S.e_q(idx,1);        % N×1
    all_e_q2{k}  = S.e_q(idx,2);        % N×1
     all_e_dq1{k}  = S.e_dq(idx,1);        % N×1
    all_e_dq2{k}  = S.e_dq(idx,2);        % N×1
    
    all_tau1{k}  = S.tau_mat(idx,1);    % N×1
    all_tau2{k}  = S.tau_mat(idx,2);    % N×1

     all_u1{k}  = S.u_mat(idx,1);    % N×1
    all_u2{k}  = S.u_mat(idx,2);    % N×1

    all_q1{k}    = S.q_use(idx,1);      % N×1
    all_q2{k}    = S.q_use(idx,2);      % N×1

    all_qd1{k}   = S.qd_mat(idx,1);     % N×1
    all_qd2{k}   = S.qd_mat(idx,2);     % N×1

     all_dq1{k}    = S.dq_use(idx,1);      % N×1
    all_dq2{k}    = S.dq_use(idx,2);      % N×1

    all_dqd1{k}   = S.dqd_mat(idx,1);     % N×1
    all_dqd2{k}   = S.dqd_mat(idx,2);     % N×1
   all_rho1{k} = S.rho_mat(idx,1);  
      all_rho2{k} = S.rho_mat(idx,2); 
%     all_rho2{k} = S.rho2(idx)';

end
labels1 = { ...
    '$x_1(0) = [\pi,\pi]$', ...
    '$x_1(0) = [\pi/2,\;\pi/2]$', ...
    '$x_1(0) = [0,\;\pi/2]$', ...
      '$x_1(0) = [\pi/4,\;\pi/4]$', ...
    '$x_1(0) = [-\pi/2,\;-\pi/2]$',  '$x_1(0) = [-3\pi/4,\;-3\pi/4]$',  'Error boundaries'  ...
};
labels2 = { ...
    '$x_1(0) = [\pi,\pi]$', ...
    '$x_1(0) = [\pi/2,\;\pi/2]$', ...
    '$x_1(0) = [0,\;\pi/2]$', ...
      '$x_1(0) = [\pi/4,\;\pi/4]$', ...
    '$x_1(0) = [-\pi/2,\;-\pi/2]$',  '$x_1(0) = [-3\pi/4,\;-3\pi/4]$',  'Desired trajectory '  ...
};

%% 画：位置跟踪误差 1
figure('Position', [100 100 1000 250]); % 设置整个figure的大小
subplot(1,2,1);
hold on;grid off;box on;
h_err = gobjects(num_cases,1);
    t = all_tspan{1}(:);
    rho1 = all_rho1{1}(:,1);  % 提取列向量
    x_fill = [t; flipud(t)];
    y_fill = [rho1; flipud(-rho1)];
fill(x_fill, y_fill,colors(1,:), ...  % 灰色
    'FaceAlpha', 0.1, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
plot(t,  rho1,  '--', 'Color',colors(1,:), 'LineWidth', 1.5)
plot(t, -rho1,  '--', 'Color', colors(1,:), 'LineWidth', 1.5)
for k = 1:num_cases
     h_err(k)=plot(all_tspan{k}, all_e_q1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(2,:));
end

    xline(3, 'LineStyle', line_styles{3},'Color',colors(4,:), 'LineWidth', 1.5, ...
        'Label', 'T=3s', 'LabelOrientation', 'horizontal', ...
        'Interpreter', 'latex');
% rho1_mat = all_rho1{1};  % N×2
% h_bound=plot(all_tspan{1}, rho1_mat(:,1), 'LineWidth',1.5,'LineStyle', line_styles{6},'Color', colors(6,:));
% plot(all_tspan{1}, rho1_mat(:,2),  'LineWidth',1.5,'LineStyle', line_styles{6},'Color', colors(6,:));
% h_all = [h_err; h_bound];
ylim([-6,1])
xlabel('$t$(s)', 'Interpreter','latex');
ylabel('$e_{1,1}$ (rad)', 'Interpreter','latex');
title('');
legend( '$\rho_1$','$-\rho_1$','$e_{1,2}$','Location', 'southeast', 'Interpreter','latex');
set(gca, 'FontSize', 12);

axMain = gca;
pos = axMain.Position;
% 插图
ax1 = axes('Position',[pos(1)+0.2*pos(3), pos(2)+0.42*pos(4), 0.3*pos(3), 0.40*pos(4)]);
box(ax1,'on'); hold(ax1,'on'); grid(ax1,'off');
plot(ax1, t,  rho1,'--','Color',colors(1,:),'LineWidth',1.5);
plot(ax1, t, -rho1,'--','Color',colors(1,:),'LineWidth',1.5);
for k=1:num_cases
    plot(ax1, all_tspan{k}, all_e_q1{k}, 'LineWidth',1.5,'LineStyle',line_styles{k},'Color',colors(2,:));
end
xlim(ax1,[0 2]); ylim(ax1,[-4.5 1.2]);
set(ax1,'FontSize',8); xlabel(ax1,''); ylabel(ax1,'');

ax2 = axes('Position',[pos(1)+0.6*pos(3), pos(2)+0.42*pos(4), 0.30*pos(3), 0.4*pos(4)]);
box(ax2,'on'); hold(ax2,'on'); grid(ax2,'off');
plot(ax2, t,  rho1,'--','Color',colors(1,:),'LineWidth',1.5);
plot(ax2, t, -rho1,'--','Color',colors(1,:),'LineWidth',1.5);
for k=1:num_cases
    plot(ax2, all_tspan{k}, all_e_q1{k}, 'LineWidth',1.5,'LineStyle',line_styles{k},'Color',colors(2,:));
end
xlim(ax2,[4 16]); ylim(ax2,[-0.25 0.25]);
set(ax2,'FontSize',8); xlabel(ax2,''); ylabel(ax2,'');
set(ax2,'Layer','top');
% 
% axes('Position',[0.2 0.35 0.1 0.4]); 
% box on;hold on;grid off;
% plot(t,  rho1,  '--', 'Color',colors(1,:), 'LineWidth', 1.5)
% plot(t, -rho1,  '--', 'Color', colors(1,:), 'LineWidth', 1.5)
% for k = 1:num_cases
%      h_err(k)=plot(all_tspan{k}, all_e_q1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(2,:));
% end
% xlim([0 2]); ylim([-4.5 1.2]);
% title('');
% 
% % 插图
% axes('Position',[0.4 0.25 0.3 0.3]); 
% box on;hold on;grid off;
% plot(t,  rho1,  '--', 'Color',colors(1,:), 'LineWidth', 1.5)
% plot(t, -rho1,  '--', 'Color', colors(1,:), 'LineWidth', 1.5)
% for k = 1:num_cases
%      h_err(k)=plot(all_tspan{k}, all_e_q1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(2,:));
% end
% xlim([4 15]); ylim([-4.5 1.2]);
% title('');







%% 画：位置跟踪误差2
subplot(1,2,2);
hold on;
h_err = gobjects(num_cases,1);
    t = all_tspan{1}(:);
    rho2 = all_rho2{1}(:,1);  % 提取列向量
    x_fill = [t; flipud(t)];
    y_fill = [rho2; flipud(-rho2)];
fill(x_fill, y_fill, colors(1,:), ...  % 灰色
    'FaceAlpha', 0.1, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
plot(t,  rho2,  '--', 'Color', colors(1,:), 'LineWidth', 1.5)
plot(t, -rho2,  '--', 'Color', colors(1,:), 'LineWidth', 1.5)

for k = 1:num_cases
     h_err(k)=plot(all_tspan{k}, all_e_q2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(2,:));
end
    xline(3, 'LineStyle', line_styles{3},'Color',colors(4,:), 'LineWidth', 1.5, ...
        'Label', 'T=3s', 'LabelOrientation', 'horizontal', ...
        'Interpreter', 'latex');
% rho1_mat = all_rho1{1};  % N×2
% h_bound=plot(all_tspan{1}, rho1_mat(:,1), '--', 'LineWidth',1.5,'LineStyle', line_styles{6},'Color', colors(6,:));
% plot(all_tspan{1}, rho1_mat(:,2), '--', 'LineWidth',1.5,'LineStyle', line_styles{6},'Color', colors(6,:));
% h_all = [h_err; h_bound];
ylim([-6,1])
xlabel('$t$(s)', 'Interpreter','latex');
ylabel('$e_{1,2}$ (rad)', 'Interpreter','latex');
title('');
% legend(h_all,labels1,'Location', 'northeast', 'Interpreter','latex');
grid off;
box on;
legend( '$\rho_2$','$-\rho_2$','$e_{1,2}$','Location', 'southeast', 'Interpreter','latex');
set(gca, 'FontSize', 12);
axMain = gca;
pos = axMain.Position;
% 插图
ax3 = axes('Position',[pos(1)+0.2*pos(3), pos(2)+0.42*pos(4), 0.3*pos(3), 0.40*pos(4)]);
box(ax3,'on'); hold(ax3,'on'); grid(ax3,'off');
plot(ax3, t,  rho2,'--','Color',colors(1,:),'LineWidth',1.5);
plot(ax3, t, -rho2,'--','Color',colors(1,:),'LineWidth',1.5);
for k=1:num_cases
    plot(ax3, all_tspan{k}, all_e_q2{k}, 'LineWidth',1.5,'LineStyle',line_styles{k},'Color',colors(2,:));
end
xlim(ax3,[0 2]); ylim(ax3,[-4.5 1.2]);
set(ax3,'FontSize',8); xlabel(ax3,''); ylabel(ax3,'');

ax4 = axes('Position',[pos(1)+0.6*pos(3), pos(2)+0.42*pos(4), 0.30*pos(3), 0.4*pos(4)]);
box(ax4,'on'); hold(ax4,'on'); grid(ax4,'off');
plot(ax4, t,  rho2,'--','Color',colors(1,:),'LineWidth',1.5);
plot(ax4, t, -rho2,'--','Color',colors(1,:),'LineWidth',1.5);
for k=1:num_cases
    plot(ax4, all_tspan{k}, all_e_q2{k}, 'LineWidth',1.5,'LineStyle',line_styles{k},'Color',colors(2,:));
end
xlim(ax4,[4 16]); ylim(ax4,[-0.2 0.2]);
xlabel(ax4,''); ylabel(ax4,'');
set(ax4,'FontSize',8); 
set(ax4,'Layer','top');

% === 导出 EPS 图像 ===
set(gcf, 'Units', 'inches', 'Position', [1 1 10 3]);

exportgraphics(gcf, 'fig5.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','none', ...
    'Resolution',600);




%% 画：关节 1 控制扭矩 tau1
figure('Position', [100 100 1000 250]); % 设置整个figure的大小
subplot(2,1,1);
hold on;
for k = 1:num_cases
     plot(all_tspan{k}, all_u1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{1},'Color',colors(1,:));
     hold on;
    plot(all_tspan{k}, all_tau1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{2},'Color',colors(2,:));
end
ylim([-32,70])
xlabel('$t$(s)', 'Interpreter','latex');
ylabel('$\tau_1$ (Nm)', 'Interpreter','latex');
legend('$\tau_1^\prime$','$\tau_1$', 'Location', 'northeast', 'Interpreter','latex');
title('');
set(gca, 'FontSize', 12);
grid off;box on;
ax5 = axes('Position',[pos(1)-0.2*pos(3), pos(2)+0.80*pos(4), 0.50*pos(3), 0.18*pos(4)]);
box(ax5,'on'); hold(ax5,'on'); grid(ax5,'off');
for k = 1:num_cases
    plot(all_tspan{k},all_u1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{1},'Color',colors(1,:));
    hold on;
    plot(all_tspan{k},  all_tau1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{2},'Color',colors(2,:));
end
xlim(ax5,[8.2 9]); ylim(ax5,[-30,30]);
xlabel(ax5,''); ylabel(ax5,'');
set(ax5,'FontSize',8); 
set(ax5,'Layer','top');
set(gcf, 'Units', 'inches', 'Position', [1 1 8 5]);
exportgraphics(gcf, 'fig6.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','none', ...
    'Resolution',600);


%% 画：关节 2 控制扭矩 tau1
subplot(2,1,2);
hold on;
for k = 1:num_cases
    plot(all_tspan{k},all_u2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{1},'Color',colors(1,:));
    hold on;
    plot(all_tspan{k},  all_tau2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{2},'Color',colors(2,:));
end
ylim([-12,28])
xlabel('$t$(s)', 'Interpreter','latex');
ylabel('$\tau_2$ (Nm)', 'Interpreter','latex');
title('');
% legend(labels1,'Location', 'northeast', 'Interpreter','latex');
grid off; box on;
legend('$\tau_2^\prime$','$\tau_2$', 'Location', 'northeast', 'Interpreter','latex');
% legend(labels1,'Location', 'northeast', 'Interpreter','latex');
title('', 'Interpreter','latex');

set(gca, 'FontSize', 12);


% 插图

ax6 = axes('Position',[pos(1)-0.2*pos(3), pos(2)+0.22*pos(4), 0.50*pos(3), 0.18*pos(4)]);
box(ax6,'on'); hold(ax6,'on'); grid(ax6,'off');
for k = 1:num_cases
    plot(all_tspan{k},all_u2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{1},'Color',colors(1,:));
    hold on;
    plot(all_tspan{k},  all_tau2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{2},'Color',colors(2,:));
end
xlim(ax6,[8.2 9]); ylim(ax6,[-12,10]);
xlabel(ax6,''); ylabel(ax6,'');
set(ax6,'FontSize',8); 
set(ax6,'Layer','top');
set(gcf, 'Units', 'inches', 'Position', [1 1 8 5]);
exportgraphics(gcf, 'fig6.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','none', ...
    'Resolution',600);


% 
