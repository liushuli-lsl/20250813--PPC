%% minimal_demo_gppf.m
clear; clc; close all;

% ---------- 1) 参数 ----------
% cfg.id = 1;
% cfg.Tp = 3;              % 预定义时间
% cfg.p  = 0.3;
% cfg.a  = 0.01;
% 
% % 后段: Σ 动力学与触发阈值（可按需微调）
% cfg.lambda_S   = 6;      % 追踪 Σ_req 速度
% cfg.k_r        = 5;      % 解除后回收强度
% cfg.Sigma_max  = 0.5;    % 最大放宽 50%
% cfg.u_d        = 0.5;    % 饱和强度归一化尺度
% cfg.theta_on   = 0.02;   % on阈
% cfg.theta_off  = 0.01;   % off阈
% cfg.delta_margin = 0.05; % 距界裕度
% cfg.gamma_u    = 1.2;    % 饱和→Σ权重
% cfg.gamma_e    = 1.0;    % 误差占用→Σ权重
% cfg.tau_u      = 0.10;   % ru LPF
% cfg.tau_e      = 0.10;   % re LPF
% cfg.rmax_Sigma = 5;      % Σ 速率限幅 (/s)
% 
% % 前段: σ（按你原逻辑）
% cfg.sigma0     = 1.0;
% cfg.k_u        = 2.0;
% cfg.k_e        = 1.0;
% cfg.sigma_min  = 0.2;
% cfg.sigma_max  = 5.0;
% cfg.ks_tight   = 1.2;
% cfg.ks_relax   = 1.8;
% cfg.rmax_sigma = 5.0;
% 配置不同通道的 G-PPF 参数（可按需微调）
cfg = struct('id',1,'Tp',T_p,'p',p,'a',a, ...   % 位置误差通道下界 a1
    'sigma0',1.0,'sigma_min',1.2,'sigma_max',4, ...
    'iota',3.0,'Sigma_max',0.5, ...
    'k_u',1.5,'k_d',2.5,'k_e',0.4, ...
    'use_lpf',true,'tau_u',0.2,'tau_d',0.2,'tau_e',0.4);
% ---------- 2) 仿真网格 ----------
dt = 1e-3;  T = 30;    % 30s
tspan = 0:dt:T;
N = numel(tspan);
n = 1;                 % 单通道演示

% ---------- 3) 构造误差和饱和残差信号（可换成你的真实数据） ----------
e      = zeros(N,1);
DeltaU = zeros(N,1);

for k = 1:N
    t = tspan(k);

    % 基本误差: 指数衰减 + 小振荡
    e(k) = 0.6*exp(-0.8*t) + 0.02*sin(6*t);

    % 在 11-13s 逼近边界（模拟“贴管”）
    if t>=11 && t<=13
        e(k) = e(k) + 0.06*sin(8*t);
    end

    % 饱和残差：在 6-8s 和 16-20s 施加成段“饱和”
    if (t>=6 && t<=8) || (t>=16 && t<=20)
        DeltaU(k) = 1.2 + 0.3*sin(40*t);   % 峰值>u_d，确保触发
    else
        DeltaU(k) = 0.0;
    end
end

% ---------- 4) 调用 gppf 并记录 ----------
rho  = zeros(N,1);
dror = zeros(N,1);
Sigma = zeros(N,1);
Sigma_req = nan(N,1);    % 仅在 t>=Tp 才有，其他置 NaN

for k = 1:N
    t = tspan(k);
    [rho(k), dror(k), aux] = gppf(t, n, e(k), DeltaU(k), DeltaU(k),cfg);
    Sigma(k) = aux.Sigma(1);
    if isfield(aux,'Sigma_req')
        Sigma_req(k) = aux.Sigma_req(1);
    end
end

% 为观察，计算 |e|、|DeltaU|
abs_e  = abs(e);
abs_du = abs(DeltaU);

% ---------- 5) 绘图 ----------
figure('Color','w','Position',[120 80 980 520]);

subplot(3,1,1);
plot(tspan, rho, 'LineWidth',1.6); hold on;
plot(tspan, abs_e, '--', 'LineWidth',1.2);
xline(cfg.Tp,':','Tp','LabelOrientation','horizontal');
ylabel('\rho(t), |e(t)|');
ylim([-0.5,2]);
legend('\rho(t)','|e(t)|','Location','northeast'); grid on;

subplot(3,1,2);
plot(tspan, abs_du, 'LineWidth',1.6); hold on;
ylabel('|\Delta U(t)|');
legend('saturation residue','Location','northeast'); grid on;

subplot(3,1,3);
plot(tspan, Sigma, 'LineWidth',1.6); hold on;
plot(tspan, Sigma_req, '--', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('\Sigma,\ \Sigma_{req}');
legend('\Sigma(t)','\Sigma_{req}(t)','Location','northeast'); grid on;

sgtitle('GPPF demo: \rho, |e|, |\Delta U|, \Sigma, \Sigma_{req}');
