%% compare_x1_to_x5.m
clearvars; close all; clc;

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
     fname = sprintf('method%d.mat', 2);
    S = load(fname, ...
        'tspan','e_q','e_dq','tau_mat','u_mat','q_use','qd_mat','dq_use','dqd_mat','rho_mat');
    t_full = S.tspan(1:end-1);
idx    = t_full <=9;        % 逻辑索引
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

T_cutoff = 3;          % 只在 3 s 之后启用阈值
th = 2; 
%% 画：位置跟踪误差 1
figure('Position', [100 100 1000 250]); % 设置整个figure的大小
subplot(1,2,1);
hold on;
h_err = gobjects(num_cases,1);

    t = all_tspan{1}(:);
    rho1 = all_rho1{1}(:,1);  % 提取列向量
    x_fill = [t; flipud(t)];
    y_fill = [rho1; flipud(-rho1)];
fill(x_fill, y_fill,colors(1,:), ...  % 灰色
    'FaceAlpha', 0.1, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
plot(t,  rho1,  '--', 'Color',colors(1,:), 'LineWidth', 1)
plot(t, -rho1,  '--', 'Color', colors(1,:), 'LineWidth', 1)
% % 绘制四条竖线

%     xline(3, 'LineStyle', line_styles{1},'Color',colors(7,:), 'LineWidth', 1.2, ...
%         'Label', 'T=3s', 'LabelOrientation', 'horizontal', ...
%         'Interpreter', 'latex');

for k = 1:num_cases

    e_k = all_e_q1{k}(:);

    % —— 选择阈值逻辑：>0.01（若想 |e|>0.01，改成 abs(e_k)>0.01）
    idx = find(t>= T_cutoff & abs(e_k) > th, 1, 'first');
    if ~isempty(idx)
        e_plot = e_k;
        e_plot(idx+1:end) = NaN;  % 让超过阈值后的段不再绘制
        % 也可用截断：e_plot = e_k(1:idx); t_k = t_k(1:idx);
    else
        e_plot = e_k;             % 整段都不超过阈值，完整绘制
    end

     h_err(k)=plot(t, e_plot, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(2,:));
end
% rho1_mat = all_rho1{1};  % N×2
% h_bound=plot(all_tspan{1}, rho1_mat(:,1), 'LineWidth',1.5,'LineStyle', line_styles{6},'Color', colors(6,:));
% plot(all_tspan{1}, rho1_mat(:,2),  'LineWidth',1.5,'LineStyle', line_styles{6},'Color', colors(6,:));
% h_all = [h_err; h_bound];
ylim([-4.5,1])
xlabel('Time (s)');
ylabel('$z_{1}$ (rad)', 'Interpreter','latex');
title('');
legend( '$\rho_1$','$-\rho_1$','$z_1$','Location', 'southeast', 'Interpreter','latex'); 
grid off;
box on;

% 插图
ax1 = gca; 
pos = get(ax1,'Position');
% ax2 = axes('Position',[0.5 0.75 0.35 0.15]); 
ax2 = axes('Position', [pos(1)+0.55*pos(3), pos(2)+0.44*pos(4), 0.25*pos(3), 0.35*pos(4)], ...   左下宽高
           'Color','none', 'Box','on');
uistack(ax2,'top'); hold(ax2,'on');
hold(ax2,'on');
plot(t,  rho1,  '--', 'Color',colors(1,:), 'LineWidth', 1)
plot(t, -rho1,  '--', 'Color', colors(1,:), 'LineWidth', 1)
for k = 1:num_cases
    plot(ax2, all_tspan{k}, all_e_q1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(2,:));
end

% yline(ax2,  0.01,'--','Color',colors(6,:), 'LineWidth', 1.2); yline(ax2, -0.01,'--','Color',colors(6,:), 'LineWidth', 1.2);
xlim(ax2,[8.5,9]); ylim(ax2,[-0.08,0.02]);
box(ax2,'on'); grid(ax2,'off');
title(ax2,'');
grid off;box on;






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
plot(t,  rho2,  '--', 'Color', colors(1,:), 'LineWidth', 1)
plot(t, -rho2,  '--', 'Color', colors(1,:), 'LineWidth', 1)

% % 绘制四条竖线

%     xline(3, 'LineStyle', line_styles{1},'Color',colors(7,:), 'LineWidth', 1.2, ...
%         'Label', 'T=3s', 'LabelOrientation', 'horizontal', ...
%         'Interpreter', 'latex');
for k = 1:num_cases
      e_k = all_e_q2{k}(:);

    % —— 选择阈值逻辑：>0.01（若想 |e|>0.01，改成 abs(e_k)>0.01）
    idx = find(t>= T_cutoff & abs(e_k) > th, 1, 'first');
    if ~isempty(idx)
        e_plot = e_k;
        e_plot(idx+1:end) = NaN;  % 让超过阈值后的段不再绘制
        % 也可用截断：e_plot = e_k(1:idx); t_k = t_k(1:idx);
    else
        e_plot = e_k;             % 整段都不超过阈值，完整绘制
    end

     h_err(k)=plot(t, e_plot, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(2,:));
end
% rho1_mat = all_rho1{1};  % N×2
% h_bound=plot(all_tspan{1}, rho1_mat(:,1), '--', 'LineWidth',1.5,'LineStyle', line_styles{6},'Color', colors(6,:));
% plot(all_tspan{1}, rho1_mat(:,2), '--', 'LineWidth',1.5,'LineStyle', line_styles{6},'Color', colors(6,:));
% h_all = [h_err; h_bound];
ylim([-4.5,1])
xlabel('Time (s)');
ylabel('$z_{2}$ (rad)', 'Interpreter','latex');
title('');
% legend(h_all,labels1,'Location', 'northeast', 'Interpreter','latex');
grid off;
box on;
legend( '$\rho_2$','$-\rho_2$','$z_2$','Location', 'southeast', 'Interpreter','latex');

% 插图
ax1 = gca; 
pos = get(ax1,'Position');
% ax2 = axes('Position',[0.5 0.75 0.35 0.15]); 
ax2 = axes('Position', [pos(1)+0.55*pos(3), pos(2)+0.44*pos(4), 0.25*pos(3), 0.35*pos(4)], ...   左下宽高
           'Color','none', 'Box','on');
uistack(ax2,'top'); hold(ax2,'on');
hold(ax2,'on');
plot(t,  rho2,  '--', 'Color',colors(1,:), 'LineWidth', 1)
plot(t, -rho2,  '--', 'Color', colors(1,:), 'LineWidth', 1)
for k = 1:num_cases
    plot(ax2, all_tspan{k}, all_e_q2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{k},'Color',colors(2,:));
end

% yline(ax2,  0.01,'--','Color',colors(6,:), 'LineWidth', 1.2); yline(ax2, -0.01,'--','Color',colors(6,:), 'LineWidth', 1.2);
xlim(ax2,[8.5,9]); ylim(ax2,[-0.08,0.15]);
box(ax2,'on'); grid(ax2,'off');
title(ax2,'');
grid off;box on;
% === 导出 EPS 图像 ===
set(gcf, 'Units', 'inches', 'Position', [1 1 10 2.5]);
set(gca, 'FontName', 'Times New Roman');
exportgraphics(gcf, 'fig9.eps', ...
    'ContentType','vector', ...
    'BackgroundColor','none', ...
    'Resolution',600);




%% 画：关节 1 位置 vs 参考 q1 vs qd1
figure('Position', [100 100 1000 250]); % 设置整个figure的大小
subplot(1,2,1);
hold on;
for k = 1:num_cases
  plot(all_tspan{k}, all_q1{k},  'LineWidth',1.5, 'LineStyle', line_styles{1},'Color',colors(1,:));
end
plot(all_tspan{1}, all_qd1{1},  'LineWidth',1.2, 'LineStyle', line_styles{2},'Color',colors(2,:));
xlabel('Time (s)');
ylabel('$q_1$ (rad)','Interpreter','latex');
title('');
legend("Actual trajectory","Desired trajectory", 'Location', 'northeast', 'Interpreter','latex');
grid off;box on;

%% 画：关节 2 位置 vs 参考 q1 vs qd1
subplot(1,2,2);
hold on;
for k = 1:num_cases
  plot(all_tspan{k}, all_q2{k},  'LineWidth',1.5, 'LineStyle', line_styles{1},'Color',colors(1,:));
end
plot(all_tspan{1}, all_qd2{1},  'LineWidth',1.2, 'LineStyle', line_styles{2},'Color',colors(2,:));
xlabel('Time (s)');
ylabel('$q_2$ (rad)','Interpreter','latex');
title('');
% legend(labels2, 'Location', 'northeast', 'Interpreter','latex');
legend("Actual trajectory","Desired trajectory", 'Location', 'northeast', 'Interpreter','latex');
grid off;box on;

% % === 导出 EPS 图像 ===
% set(gcf, 'Units', 'inches', 'Position', [1 1 10 2.5]);
% set(gca, 'FontName', 'Times New Roman');
% exportgraphics(gcf, 'fig9.eps', ...
%     'ContentType','vector', ...
%     'BackgroundColor','none', ...
%     'Resolution',600);



%% 画：关节 1 控制扭矩 tau1
figure('Position', [100 100 1000 250]); % 设置整个figure的大小
subplot(1,2,1);
hold on;
for k = 1:num_cases
    plot(all_tspan{k}, all_u1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{1},'Color',colors(2,:));
   hold  on ;
     plot(all_tspan{k}, all_tau1{k}, 'LineWidth',1.5, 'LineStyle', line_styles{2},'Color',colors(1,:));
end
ylim([-20,32])
xlabel('Time (s)');
ylabel('$u_1$ (Nm)', 'Interpreter','latex');
legend('$u_1$', '$u_1^\prime$','Location', 'northeast', 'Interpreter','latex');
title('');

grid off;box on;


%% 画：关节 2 控制扭矩 tau1
subplot(1,2,2);
hold on;
for k = 1:num_cases
    
     
    plot(all_tspan{k}, all_u2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{1},'Color',colors(2,:));
      hold on;
    plot(all_tspan{k}, all_tau2{k}, 'LineWidth',1.5, 'LineStyle', line_styles{2},'Color',colors(1,:));
end
ylim([-10,12])
xlabel('Time (s)');
ylabel('$u_2$ (Nm)', 'Interpreter','latex');
title('');
% legend(labels1,'Location', 'northeast', 'Interpreter','latex');
grid off; box on;
legend('$u_2$', '$u_2^\prime$','Location', 'northeast', 'Interpreter','latex');
% legend(labels1,'Location', 'northeast', 'Interpreter','latex');
title('', 'Interpreter','latex');

% % === 导出 EPS 图像 ===
% set(gcf, 'Units', 'inches', 'Position', [1 1 10 2.5]);
% set(gca, 'FontName', 'Times New Roman');
% exportgraphics(gcf, 'fig3.eps', ...
%     'ContentType','vector', ...
%     'BackgroundColor','none', ...
%     'Resolution',600);


function mark_point(ax, x0, y0, varargin)
%MARK_POINT  在坐标轴 ax 上按数据坐标标注圆圈+箭头+文本
% mark_point(ax, x0, y0, 'Radius', r, 'ArrowFrom', [xa ya], 'Text', '说明', 'Color', 'r')
%
% 必要参数：
%   ax        目标坐标轴句柄（如 gca 或 subplot 返回值）
%   x0,y0     需要圈出的点（数据坐标）
%
% 可选参数（Name-Value）：
%   'Radius'    圆半径，默认根据轴范围取 2%（数据尺度）
%   'ArrowFrom' 箭头起点 [xa ya]（数据坐标）。若省略则自动放在右上偏移处
%   'Text'      文本内容，默认 'key point'
%   'Color'     颜色，默认 'r'
%   'LineWidth' 线宽，默认 1.5
%   'TextOffset' 文本相对箭头起点的微调 [dx dy]（数据坐标），默认 [0 0]

% 解析参数
p = inputParser;
addParameter(p, 'Radius', [], @(v) isempty(v) || isscalar(v));
addParameter(p, 'ArrowFrom', [], @(v) isempty(v) || (isnumeric(v)&&numel(v)==2));
addParameter(p, 'Text', 'key point', @(s) isstring(s)||ischar(s));
addParameter(p, 'Color', 'r');
addParameter(p, 'LineWidth', 1.5, @isscalar);
addParameter(p, 'TextOffset', [0 0], @(v) isnumeric(v)&&numel(v)==2);
parse(p, varargin{:});
R   = p.Results.Radius;
P0  = [x0 y0];
Pa  = p.Results.ArrowFrom;
txt = p.Results.Text;
col = p.Results.Color;
lw  = p.Results.LineWidth;
toff= p.Results.TextOffset;

hold(ax,'on');

% 自动半径（按轴范围的 2%）
if isempty(R)
    xr = diff(xlim(ax)); yr = diff(ylim(ax));
    R = 0.02 * max(xr, yr);
end

% 1) 画红色空心圆（数据坐标）
rectangle(ax, 'Position',[x0-R, y0-R, 2*R, 2*R], ...
    'Curvature',[1 1], 'EdgeColor',col, 'LineWidth',lw, 'Clipping','off');

% 2) 箭头起点（若未给，则放在右上偏移处）
if isempty(Pa)
    xr = diff(xlim(ax)); yr = diff(ylim(ax));
    Pa = P0 + [0.12*xr, 0.12*yr];   % 右上角偏移
end

% 3) 用 quiver 在数据坐标画箭头
dvec = P0 - Pa;  % 指向 (x0,y0)
quiver(ax, Pa(1), Pa(2), dvec(1), dvec(2), 0, ...
    'Color',col,'LineWidth',lw,'MaxHeadSize',0.5);

% 4) 文本（放在箭头起点附近）
text(ax, Pa(1)+toff(1), Pa(2)+toff(2), ['  ' txt], ...
    'Color',col,'FontWeight','bold','Interpreter','latex', ...
    'VerticalAlignment','middle','Clipping','off');
end

