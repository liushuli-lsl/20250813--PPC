%% compare_x1_to_x5.m
clearvars; close all; clc;

num_cases = 4;
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

for k = 1:num_cases
    fname = sprintf('T%d.mat', k);
    S = load(fname, ...
        'tspan','e_q','e_dq','tau_mat','q_use','qd_mat','dq_use','dqd_mat','rho1','rho2');
    t_full = S.tspan(1:end-1);
idx    = t_full <= 6;        % 逻辑索引
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
        r1 = S.rho1(idx);   % N×1
    r2 = S.rho2(idx);   % N×1

    % 扩展成 N×2，每列都是同一个性能函数
    all_rho1{k} = [r1, r1];  
    all_rho2{k} = [r2, r2];

end
labels1 = { ...
  
    '$Tp = 1$s', ...
    '$Tp = 2$s', ...
    '$Tp = 3$s', ...
    '$Tp = 5$s', ...
%      'Predifined Time'  ...
};
labels2 = { ...
  
    '$Tp = 1$s', ...
    '$Tp = 2$s', ...
    '$Tp= 3$s', ...
    '$Tp = 5$s', 'Desired trajectory'  ...
};
line_positions = [1, 3, 5]; % 四条竖线的x坐标位置
%% 画：位置跟踪误差 1
figure('Position', [100 100 800 600]); % 设置整个figure的大小
subplot(2,1,1);
hold on;
h_err = gobjects(num_cases,1);
for k = 1:num_cases
     h_err(k)=plot(all_tspan{k}, all_e_q1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end
% rho1_mat = all_rho1{1};  % N×2
% h_bound=plot(all_tspan{1}, rho1_mat(:,1), 'LineStyle', line_styles{5},'Color', colors(6,:));
% plot(all_tspan{1}, -rho1_mat(:,2),  'LineWidth',1.5,'LineStyle', line_styles{5},'Color', colors(6,:));
% h_all = [h_err; h_bound];

% % 绘制四条竖线
% for i = 1:length(line_positions)
%     xline(line_positions(i), '--', 'LineStyle', line_styles{1},'Color',colors(5,:), 'LineWidth', 1.2, ...
%         'Label', labels1{i}, 'LabelOrientation', 'horizontal', ...
%         'Interpreter', 'latex');
% end
yline(0.01,  'Label', '0.01', 'LabelHorizontalAlignment', 'right',  'LineStyle', line_styles{2},'Color',colors(5,:), 'LineWidth', 1.2);  % 标签靠右
yline(-0.01, 'Label', '-0.01', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LineStyle', line_styles{2},'Color',colors(5,:), 'LineWidth', 1.2);  % 标签靠右

ylim([-0.25,0.1])
xlabel('Time (s)');
ylabel('Position error e_{1} (rad)');
title('');
legend(labels1,'Location', 'southeast', 'Interpreter','latex');
grid off;

%% 画：位置跟踪误差2
subplot(2,1,2);
hold on;
h_err = gobjects(num_cases,1);
for k = 1:num_cases
     h_err(k)=plot(all_tspan{k}, all_e_q2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end
% rho1_mat = all_rho1{1};  % N×2
% h_bound=plot(all_tspan{1}, rho2_mat(:,1), '--', 'LineStyle', line_styles{k},'Color', colors(6,:));
% plot(all_tspan{1}, -rho2_mat(:,2), '--', 'LineWidth',1.5,'LineStyle', line_styles{k},'Color', colors(6,:));
% h_all = [h_err; h_bound];
% 绘制四条竖线
% for i = 1:length(line_positions)
%     xline(line_positions(i), '--', 'LineStyle', line_styles{1},'Color',colors(5,:), 'LineWidth', 1.2, ...
%         'Label', labels1{i}, 'LabelOrientation', 'horizontal', ...
%         'Interpreter', 'latex');
% end
yline(0.01,  'Label', '0.01', 'LabelHorizontalAlignment', 'right',  'LineStyle', line_styles{2},'Color',colors(5,:), 'LineWidth', 1.2);  % 标签靠右
yline(-0.01, 'Label', '-0.01', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LineStyle', line_styles{2},'Color',colors(5,:), 'LineWidth', 1.2);  % 标签靠右
xlabel('Time (s)');
ylabel('Position error e_{2} (rad)');
title('');
legend(labels1,'Location', 'southeast', 'Interpreter','latex');
grid off;




%% 画：速度跟踪误差 e_q1
figure('Position', [100 100 800 600]); % 设置整个figure的大小
subplot(2,1,1);
hold on;
h_err = gobjects(num_cases,1);
for k = 1:num_cases
     h_err(k)=plot(all_tspan{k}, all_e_dq1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end
% rho1_mat = all_rho1{1};  % N×2
% h_bound=plot(all_tspan{1}, rho1_mat(:,1), '--', 'LineStyle', line_styles{k},'Color', colors(6,:));
% plot(all_tspan{1}, -rho1_mat(:,2), '--', 'LineWidth',1.5,'LineStyle', line_styles{k},'Color', colors(6,:));
% h_all = [h_err; h_bound];

% 绘制四条竖线
% for i = 1:length(line_positions)
%     xline(line_positions(i), '--', 'LineStyle', line_styles{1},'Color',colors(5,:), 'LineWidth', 1.2, ...
%         'Label', labels1{i}, 'LabelOrientation', 'horizontal', ...
%         'Interpreter', 'latex');
% end
yline(0.2,  'Label', '0.2', 'LabelHorizontalAlignment', 'right',  'LineStyle', line_styles{2},'Color',colors(5,:), 'LineWidth', 1.2);  % 标签靠右
yline(-0.2, 'Label', '-0.2', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LineStyle', line_styles{2},'Color',colors(5,:), 'LineWidth', 1.2);  % 标签靠右

ylim([-1,5])
xlabel('Time (s)');
ylabel('Velocity error $\dot{e}_{1}$ (rad/s)','Interpreter','latex');
title('');
legend(labels1,'Location', 'northeast', 'Interpreter','latex');
grid off;

%% 画：速度跟踪误差 e_q2
subplot(2,1,2);
hold on;
h_err = gobjects(num_cases,1);
for k = 1:num_cases
     h_err(k)=plot(all_tspan{k}, all_e_dq2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end
% rho1_mat = all_rho1{1};  % N×2
% h_bound=plot(all_tspan{1}, rho2_mat(:,1), '--', 'LineStyle', line_styles{k},'Color', colors(6,:));
% plot(all_tspan{1}, -rho2_mat(:,2), '--', 'LineWidth',1.5,'LineStyle', line_styles{k},'Color', colors(6,:));
% h_all = [h_err; h_bound];
% 绘制四条竖线
% for i = 1:length(line_positions)
%     xline(line_positions(i), '--', 'LineStyle', line_styles{1},'Color',colors(5,:), 'LineWidth', 1.2, ...
%         'Label', labels1{i}, 'LabelOrientation', 'horizontal', ...
%         'Interpreter', 'latex');
% end
yline(0.2,  'Label', '0.2', 'LabelHorizontalAlignment', 'right',  'LineStyle', line_styles{2},'Color',colors(5,:), 'LineWidth', 1.2);  % 标签靠右
yline(-0.2, 'Label', '-0.2', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'LineStyle', line_styles{2},'Color',colors(5,:), 'LineWidth', 1.2);  % 标签靠右
ylim([-2,2])
xlabel('Time (s)');
ylabel('Velocity error $\dot{e}_{2}$ (rad/s)','Interpreter','latex');
title('');
legend( labels1,'Location', 'northeast', 'Interpreter','latex');
grid off;





%% 画：关节 1 位置 vs 参考 q1 vs qd1
figure('Position', [100 100 800 600]); % 设置整个figure的大小
subplot(2,1,1);
hold on;
for k = 1:num_cases
  plot(all_tspan{k}, all_q1{k}, '-', 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end
plot(all_tspan{1}, all_qd1{1}, '--', 'LineWidth',1.2, 'LineStyle', line_styles{5},'Color',colors(5,:));
xlabel('Time (s)');
ylabel('Position tracking $q_1$ (rad)','Interpreter','latex');
title('');
legend(labels2, 'Location', 'northeast', 'Interpreter','latex');
grid off;

%% 画：关节 2 位置 vs 参考 q1 vs qd1
subplot(2,1,2);
hold on;
for k = 1:num_cases
  plot(all_tspan{k}, all_q2{k}, '-', 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end
plot(all_tspan{1}, all_qd2{2}, '--', 'LineWidth',1.2, 'LineStyle', line_styles{5},'Color',colors(5,:));
xlabel('Time (s)');
ylabel('Position tracking $q_2$ (rad)','Interpreter','latex');
title('');
legend(labels2, 'Location', 'northeast', 'Interpreter','latex');
grid off;





%% 画：关节 1 速度 vs 参考 q1 vs qd1
figure('Position', [100 100 800 600]); % 设置整个figure的大小
subplot(2,1,1);
hold on;
for k = 1:num_cases
  plot(all_tspan{k}, all_dq1{k}, '-', 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end
plot(all_tspan{1}, all_dqd1{1}, '--', 'LineWidth',1.2, 'LineStyle', line_styles{5},'Color',colors(5,:));
xlabel('Time (s)');
ylabel('Velocity tracking $\dot{q}_1$ (rad/s)','Interpreter','latex');
title('');

legend(labels2, 'Location', 'northeast', 'Interpreter','latex');
grid off;

%% 画：关节 2 速度 vs 参考 q1 vs qd1
subplot(2,1,2);
hold on;
for k = 1:num_cases
  plot(all_tspan{k}, all_q2{k}, '-', 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end
plot(all_tspan{1}, all_qd2{2}, '--', 'LineWidth',1.2, 'LineStyle', line_styles{5},'Color',colors(5,:));
xlabel('Time (s)');
ylabel('Velocity tracking $\dot{q}_2$ (rad/s)','Interpreter','latex');
title('');

legend(labels2, 'Location', 'northeast', 'Interpreter','latex');
grid off;






%% 画：关节 1 控制扭矩 tau1
figure('Position', [100 100 800 600]); % 设置整个figure的大小
subplot(2,1,1);
hold on;
for k = 1:num_cases
    plot(all_tspan{k}, all_tau1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end
ylim([-20,20])
xlabel('Time (s)');
ylabel('Torque u_1 (Nm)');
title('');
legend(labels1,'Location', 'northeast', 'Interpreter','latex');
grid off;


%% 画：关节 2 控制扭矩 tau1
subplot(2,1,2);
hold on;
for k = 1:num_cases
    plot(all_tspan{k}, all_tau2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(k,:));
end
ylim([-20,20])
xlabel('Time (s)');
ylabel('Torque $u_2$ (Nm)');
title('');
legend(labels1,'Location', 'northeast', 'Interpreter','latex');
grid off;


