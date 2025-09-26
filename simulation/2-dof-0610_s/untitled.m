
% 演示并绘制：rho(t), drho(t), drho/rho, 以及 ru, rd, re, g(t), Sigma(t)
% 依赖：gppf.m（使用你贴的版本）

clear; close all; clc;

%% ===== 配置（与控制器一致）=====
Tp = 3.0;  p = 0.3;
cfg = struct('id', 101, 'Tp', Tp, 'p', p, 'a', 0.02, ...
             'sigma0', 1.0, 'sigma_min', 0.6, 'sigma_max', 1.5, ...
             'iota', 3.0, 'Sigma_max', 0.40, ...
             'k_u', 0.7, 'k_d', 0.4, 'k_e', 0.2, ...
             'use_lpf', true, 'tau_u', 0.6, 'tau_d', 0.8, 'tau_e', 0.8);

% 时间轴
Tend = 6.0; N = 6001;
t = linspace(0, Tend, N)';  dt = t(2)-t(1);

%（示例）误差 e、饱和残差 Δu、扰动 d（两通道）
n = 2;
e      = zeros(N, n);
DeltaU = zeros(N, n);
dd     = zeros(N, n);

% 构造：e(t) 平滑衰减 + 微小振荡（两通道不同幅度）
e(:,1) = 0.6*exp(-2.0*t) + 0.05*sin(4*t);
e(:,2) = 0.4*exp(-1.2*t) + 0.04*sin(3*t);

% 构造：Δu 与 d 在 Tp 之后出现两个高斯“事件”，模拟饱和与扰动
gauss = @(x,mu,s) exp(-((x-mu)./s).^2);
DeltaU(:,1) = 0.5*gauss(t, 3.6, 0.25) + 0.35*gauss(t, 4.8, 0.30);
DeltaU(:,2) = 0.4*gauss(t, 3.9, 0.28) + 0.25*gauss(t, 5.1, 0.35);
dd(:,1)     = 0.35*gauss(t, 3.9, 0.35) + 0.20*gauss(t, 5.2, 0.35);
dd(:,2)     = 0.30*gauss(t, 4.0, 0.30) + 0.18*gauss(t, 5.0, 0.40);

%% ===== 预分配日志 =====
rho   = zeros(N, n);
drho  = zeros(N, n);
ratio = zeros(N, n);
ru    = zeros(N, n);
rd    = zeros(N, n);
re    = zeros(N, n);
gg    = zeros(N, 1);
Sigma = zeros(N, n);

%% ===== 主循环：先算 (rho,drho,ratio)，再回灌 feed_only 更新LPF =====
% 建议：第一次运行前清除函数的 persistent 状态
clear gppf

for k = 1:N
    ek   = e(k,:).';
    duk  = DeltaU(k,:).';
    dk   = dd(k,:).';

    % 计算 rho, drho, ratio；本步用当前的 LPF 状态
    [r, dr, aux] = gppf(t(k), ek, duk*0, dk, cfg);  % 先用0占位Δu（避免环）
    rho(k,:)   = r(:).';
    drho(k,:)  = dr(:).';
    ratio(k,:) = aux.ratio(:).';
    ru(k,:)    = aux.ru(:).';
    rd(k,:)    = aux.rd(:).';
    re(k,:)    = aux.re(:).';
    gg(k)      = aux.g;
    Sigma(k,:) = aux.Sigma(:).';

    % 回灌“真实” Δu 与 d：只更新 LPF，不改变当步 rho（避免代数环）
    cfg.feed_only = true;
    gppf(t(k), ek, duk, dk, cfg);
    cfg.feed_only = false;
end

%% ===== 作图 =====
lw = 1.7;
figure('Color','w','Position',[80 80 1000 720]);

subplot(2,2,1); hold on; box on;
plot(t,  rho(:,1), 'LineWidth',lw);
plot(t,  rho(:,2), 'LineWidth',lw);
yline(cfg.a,'k--','a','HandleVisibility','off');
xline(Tp,'k--','T_p','LabelVerticalAlignment','bottom'); 
title('\rho(t)'); xlabel('t (s)'); ylabel('\rho');
legend('\rho_1','\rho_2','Location','best'); grid on;

subplot(2,2,2); hold on; box on;
plot(t, drho(:,1), 'LineWidth',lw);
plot(t, drho(:,2), 'LineWidth',lw);
xline(Tp,'k--','T_p');
title('\dot{\rho}(t)'); xlabel('t (s)'); ylabel('\dot{\rho}');
legend('\dot{\rho}_1','\dot{\rho}_2','Location','best'); grid on;

subplot(2,2,3); hold on; box on;
plot(t, ratio(:,1), 'LineWidth',lw);
plot(t, ratio(:,2), 'LineWidth',lw);
xline(Tp,'k--','T_p');
title('\dot{\rho}/\rho (for BLF compensation)'); xlabel('t (s)'); ylabel('\dot{\rho}/\rho');
legend('ratio_1','ratio_2','Location','best'); grid on;

subplot(2,2,4); hold on; box on;
plot(t, ru(:,1), 'LineWidth',lw); 
plot(t, rd(:,1), 'LineWidth',lw); 
plot(t, re(:,1), 'LineWidth',lw);
plot(t, gg,      'LineWidth',lw+0.2);
plot(t, Sigma(:,1), 'LineWidth',lw+0.3);
xline(Tp,'k--','T_p');
title('Drivers (ch-1) & gate'); xlabel('t (s)');
legend('r_u','r_d','r_e','g(t)','\Sigma','Location','best'); grid on;

sgtitle('G-PPF key curves (2 channels)');


