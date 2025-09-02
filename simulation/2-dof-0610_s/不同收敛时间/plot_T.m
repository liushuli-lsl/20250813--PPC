%% compare_x1_to_x5.m
clearvars; close all; clc;

num_cases = 3;
% tp_list = [1, 3, 5];
% num_cases = numel(tp_list);
colors = lines(6);
line_styles = {'-', '--', '-.', ':','-', '--',}; % 不同线型

% 准备存储
all_tspan = cell(num_cases,1);
all_e_q1   = cell(num_cases,1);
all_e_q2   = cell(num_cases,1);
all_e_dq1   = cell(num_cases,1);
all_e_dq2   = cell(num_cases,1);
all_tau1   = cell(num_cases,1);
all_tau2   = cell(num_cases,1);
all_q1     = cell(num_cases,1);
all_q2     = cell(num_cases,1);
all_qd1    = cell(num_cases,1);
all_qd2    = cell(num_cases,1);
all_dq1     = cell(num_cases,1);
all_dq2     = cell(num_cases,1);
all_dqd1    = cell(num_cases,1);
all_dqd2    = cell(num_cases,1);
 all_rho1 = cell(num_cases,1);
    all_rho2 = cell(num_cases,1);

for k = 1:num_cases
    fname = sprintf('T%d.mat', k);
    S = load(fname, ...
        'tspan','e_q','e_dq','tau_mat','q_use','qd_mat','dq_use','dqd_mat','rho1','rho2');
    t_full = S.tspan(1:end-1);
idx    = t_full <= 5;        % 逻辑索引
    all_tspan{k} = S.tspan(idx);           % 1×N
    all_e_q1{k}  = S.e_q(idx,1);        % N×1
    all_e_q2{k}  = S.e_q(idx,2);        % N×1
     all_e_dq1{k}  = S.e_dq(idx,1);        % N×1
    all_e_dq2{k}  = S.e_dq(idx,2);        % N×1
    all_tau1{k}  = S.tau_mat(idx,1);    % N×1
    all_tau2{k}  = S.tau_mat(idx,2);    % N×1
    all_q1{k}    = S.q_use(idx,1);      % N×1
    all_q2{k}    = S.q_use(idx,2);      % N×1
    all_qd1{k}   = S.qd_mat(idx,1);     % N×1
    all_qd2{k}   = S.qd_mat(idx,2);     % N×1
     all_dq1{k}    = S.dq_use(idx,1);      % N×1
    all_dq2{k}    = S.dq_use(idx,2);      % N×1
    all_dqd1{k}   = S.dqd_mat(idx,1);     % N×1
    all_dqd2{k}   = S.dqd_mat(idx,2);     % N×1
    all_rho1{k} = S.rho1(idx)';  
    all_rho2{k} = S.rho2(idx)';
%         r1 = S.rho1(idx);   % N×1
%     r2 = S.rho2(idx);   % N×1
% 
%     % 扩展成 N×2，每列都是同一个性能函数
%     all_rho1{k} = [r1, r1];  
%     all_rho2{k} = [r2, r2];

end
labels1 = { ...
  
    '$Tp = 1.5$s', ...
    '$Tp = 2.5$s', ...
    '$Tp = 3.5$s', ...
      '$Tp = 4.5$s', ...
%      'Predifined Time'  ...
};
labels2 = { ...
  '$Tp = 1.5$s', ...
    '$Tp = 2.5$s', ...
    '$Tp = 3.5$s', ...
      '$Tp = 4.5$s', ... 'Desired trajectory'  ...
};
line_positions = [1.5, 2.5, 3.5]; % 四条竖线的x坐标位置

%% 画：位置跟踪误差 1
figure('Position', [100 100 600 600]); % 设置整个figure的大小
subplot(1,2,1);
hold on;
for k = 1:num_cases
    t = all_tspan{k}(:);
    rho1 = all_rho1{k}(:,1);  % 提取列向量
    x_fill = [t; flipud(t)];
    y_fill = [rho1; flipud(-rho1)];

fill(x_fill, y_fill, [0.5 0.5 0.5], ...  % 灰色
    'FaceAlpha', 0.06, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
plot(t,  rho1,  '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8)
plot(t, -rho1,  '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8)
end
% % 绘制四条竖线
for i = 1:length(line_positions)
    xline(line_positions(i), 'LineStyle', line_styles{1},'Color',colors(4,:), 'LineWidth', 1.2, ...
        'Label', labels1{i}, 'LabelOrientation', 'horizontal', ...
        'Interpreter', 'latex');
end
% === 再画误差曲线 ===
h_err = gobjects(num_cases,1);
for k = 1:num_cases
    h_err(k)=plot(all_tspan{k}, all_e_q1{k}, ...
        'LineWidth',1.5, ...
        'LineStyle', line_styles{k}, ...
        'Color', colors(k,:));
end

% === 图例绑定曲线（不包括填充区域）===

% yline(0.01,  'Label', '0.01', 'LabelHorizontalAlignment', 'right',  'LineStyle', line_styles{5},'Color',colors(6,:), 'LineWidth', 1.2);  % 标签靠右
% yline(-0.01, 'Label', '-0.01', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LineStyle', line_styles{5},'Color',colors(6,:), 'LineWidth', 1.2);  % 标签靠右
ylim([-0.2,0.8])
xlabel('Time (s)');
ylabel('$e_{1}$ (rad)', 'Interpreter','latex');
title('');
grid off;box on;


%% 画：位置跟踪误差2
subplot(1,2,2);
hold on;
h_err = gobjects(num_cases,1);
for k = 1:num_cases
    t = all_tspan{k}(:);
    rho1 = all_rho1{k}(:,1);  % 提取列向量
    x_fill = [t; flipud(t)];
    y_fill = [rho1; flipud(-rho1)];

fill(x_fill, y_fill, [0.5 0.5 0.5], ...  % 灰色
    'FaceAlpha', 0.06, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
plot(t,  rho1,  '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8)
plot(t, -rho1,  '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8)
end
% % 绘制四条竖线
for i = 1:length(line_positions)
    xline(line_positions(i), 'LineStyle', line_styles{1},'Color',colors(4,:), 'LineWidth', 1.2, ...
        'Label', labels1{i}, 'LabelOrientation', 'horizontal', ...
        'Interpreter', 'latex');
end
for k = 1:num_cases
     h_err(k)=plot(all_tspan{k}, all_e_q2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end

% 绘制四条竖线
% for i = 1:length(line_positions)
%     xline(line_positions(i), '--', 'LineStyle', line_styles{1},'Color',colors(6,:), 'LineWidth', 1.2, ...
%         'Label', labels1{i}, 'LabelOrientation', 'horizontal', ...
%         'Interpreter', 'latex');
% end
% yline(0.01,  'Label', '0.01', 'LabelHorizontalAlignment', 'right',  'LineStyle', line_styles{5},'Color',colors(6,:), 'LineWidth', 1.2);  % 标签靠右
% yline(-0.01, 'Label', '-0.01', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LineStyle', line_styles{5},'Color',colors(6,:), 'LineWidth', 1.2);  % 标签靠右
ylim([-0.5,0.8])
xlabel('Time (s)');
ylabel('$e_{2}$ (rad)', 'Interpreter','latex');
title('');
% legend(h_err, labels1, 'Location', 'northeast', 'Interpreter','latex');
grid off;box on;
legend(h_err, labels1, 'Location', 'southeast', 'Interpreter','latex');
set(gcf, 'Units', 'inches', 'Position', [1 1 10 2.5]);
set(gca, 'FontName', 'Times New Roman');
% exportgraphics(gcf, 'fig3.eps', ...
%     'ContentType','vector', ...
%      'BackgroundColor','white', ...
%     'Resolution',600);
exportgraphics(gcf, 'fig3.eps', 'Resolution', 1200);



%% 画：速度跟踪误差 e_q1
figure('Position', [100 100 600 600]); % 设置整个figure的大小
subplot(1,2,1);
hold on;
h_err = gobjects(num_cases,1);
for k = 1:num_cases
    t = all_tspan{k}(:);
    rho2 = all_rho2{k}(:,1);  % 提取列向量
    x_fill = [t; flipud(t)];
    y_fill = [rho2; flipud(-rho2)];

fill(x_fill, y_fill, [0.5 0.5 0.5], ...  % 灰色
    'FaceAlpha', 0.06, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
plot(t,  rho2,  '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8)
plot(t, -rho2,  '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8)
end
% % 绘制四条竖线
for i = 1:length(line_positions)
    xline(line_positions(i), 'LineStyle', line_styles{1},'Color',colors(4,:), 'LineWidth', 1.2, ...
        'Label', labels1{i}, 'LabelOrientation', 'horizontal', ...
        'Interpreter', 'latex');
end
for k = 1:num_cases
     h_err(k)=plot(all_tspan{k}, all_e_dq1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end

% yline(0.1,  'Label', '0.1', 'LabelHorizontalAlignment', 'right',  'LineStyle', line_styles{5},'Color',colors(6,:), 'LineWidth', 1.2);  % 标签靠右
% yline(-0.1, 'Label', '-0.1', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LineStyle', line_styles{5},'Color',colors(6,:), 'LineWidth', 1.2);  % 标签靠右

ylim([-3,1.5])
xlabel('Time (s)');
ylabel('$\dot{e}_{1}$ (rad/s)','Interpreter','latex');
title('');
grid off;box on;


%% 画：速度跟踪误差 e_q2
subplot(1,2,2);
hold on;
h_err = gobjects(num_cases,1);
for k = 1:num_cases
    t = all_tspan{k}(:);
    rho2 = all_rho2{k}(:,1);  % 提取列向量
    x_fill = [t; flipud(t)];
    y_fill = [rho2; flipud(-rho2)];

fill(x_fill, y_fill, [0.5 0.5 0.5], ...  % 灰色
    'FaceAlpha', 0.06, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
plot(t,  rho2,  '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8)
plot(t, -rho2,  '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8)
end
% % 绘制四条竖线
for i = 1:length(line_positions)
    xline(line_positions(i), 'LineStyle', line_styles{1},'Color',colors(4,:), 'LineWidth', 1.2, ...
        'Label', labels1{i}, 'LabelOrientation', 'horizontal', ...
        'Interpreter', 'latex');
end
for k = 1:num_cases
     h_err(k)=plot(all_tspan{k}, all_e_dq2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end

% yline(0.1,  'Label', '0.1', 'LabelHorizontalAlignment', 'right',  'LineStyle', line_styles{5},'Color',colors(6,:), 'LineWidth', 1.2);  % 标签靠右
% yline(-0.1, 'Label', '-0.1', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LineStyle', line_styles{5},'Color',colors(6,:), 'LineWidth', 1.2);  % 标签靠右
ylim([-3,1.5])
xlabel('Time (s)');
ylabel('$\dot{e}_{2}$ (rad/s)','Interpreter','latex');
title('');
% legend( labels1,'Location', 'northeast', 'Interpreter','latex');
grid off;box on;
legend(h_err, labels1, 'Location', 'southeast', 'Interpreter','latex');
set(gcf, 'Units', 'inches', 'Position', [1 1 10 2.5]);
set(gca, 'FontName', 'Times New Roman');
% exportgraphics(gcf, 'fig4.eps', ...
%     'ContentType','vector', ...
%     'BackgroundColor','none', ...
%     'Resolution',600);
exportgraphics(gcf, 'fig4.eps', 'Resolution', 1200);



% %% 画：关节 1 位置 vs 参考 q1 vs qd1
% figure('Position', [100 100 600 300]); % 设置整个figure的大小
% subplot(2,1,1);
% hold on;
% for k = 1:num_cases
%   plot(all_tspan{k}, all_q1{k}, '-', 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
% end
% plot(all_tspan{1}, all_qd1{1}, '--', 'LineWidth',1.2, 'LineStyle', line_styles{5},'Color',colors(6,:));
% xlabel('Time (s)');
% ylabel('$q_1$ (rad)','Interpreter','latex');
% title('');
% legend(labels2, 'Location', 'northeast', 'Interpreter','latex');
% grid off;box on;
% 
% %% 画：关节 2 位置 vs 参考 q1 vs qd1
% subplot(2,1,2);
% hold on;
% for k = 1:num_cases
%   plot(all_tspan{k}, all_q2{k}, '-', 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
% end
% plot(all_tspan{1}, all_qd2{2}, '--', 'LineWidth',1.2, 'LineStyle', line_styles{5},'Color',colors(6,:));
% xlabel('Time (s)');
% ylabel('$q_2$ (rad)','Interpreter','latex');
% title('');
% legend(labels2, 'Location', 'northeast', 'Interpreter','latex');
% grid off;box on;
% 
% 
% 
% 
% 
% %% 画：关节 1 速度 vs 参考 q1 vs qd1
% figure('Position', [100 100 600 300]); % 设置整个figure的大小
% subplot(2,1,1);
% hold on;
% for k = 1:num_cases
%   plot(all_tspan{k}, all_dq1{k}, '-', 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
% end
% plot(all_tspan{1}, all_dqd1{1}, '--', 'LineWidth',1.2, 'LineStyle', line_styles{5},'Color',colors(6,:));
% xlabel('Time (s)');
% ylabel('$\dot{q}_1$ (rad/s)','Interpreter','latex');
% title('');
% 
% legend(labels2, 'Location', 'northeast', 'Interpreter','latex');
% grid off;box on;
% 
% %% 画：关节 2 速度 vs 参考 q1 vs qd1
% subplot(2,1,2);
% hold on;
% for k = 1:num_cases
%   plot(all_tspan{k}, all_q2{k}, '-', 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
% end
% plot(all_tspan{1}, all_qd2{2}, '--', 'LineWidth',1.2, 'LineStyle', line_styles{5},'Color',colors(6,:));
% xlabel('Time (s)');
% ylabel('$\dot{q}_2$ (rad/s)','Interpreter','latex');
% title('');
% 
% legend(labels2, 'Location', 'northeast', 'Interpreter','latex');
% grid off;box on;






%% 画：关节 1 控制扭矩 tau1
figure('Position', [100 100 600 600]); % 设置整个figure的大小
subplot(1,2,1);
hold on;
for k = 1:num_cases
    plot(all_tspan{k}, all_tau1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end
% ylim([-12,12])
xlabel('Time (s)');
ylabel('$u_1$ (Nm)', 'Interpreter','latex');
title('');

grid off;box on;


%% 画：关节 2 控制扭矩 tau1
subplot(1,2,2);
hold on;
for k = 1:num_cases
    plot(all_tspan{k}, all_tau2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end
% ylim([-12,12])
xlabel('Time (s)');
ylabel('$u_2$ (Nm)', 'Interpreter','latex');
title('');
legend(labels1,'Location', 'southeast', 'Interpreter','latex');
% legend(labels1,'Location', 'northeast', 'Interpreter','latex');
grid off;box on;
set(gcf, 'Units', 'inches', 'Position', [1 1 10 2.5]);
set(gca, 'FontName', 'Times New Roman');
% print(gcf, '-depsc2', '-r600', 'fig5.eps');
% exportgraphics(gcf, 'fig5.eps', ...
%     'ContentType','vector', ...
%     'BackgroundColor','none', ...
%     'Resolution',600);
exportgraphics(gcf, 'fig5.eps', 'Resolution', 1200);

