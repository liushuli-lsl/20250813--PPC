%% === 普通性能函数 PPF vs 全局自适应性能函数 G-PPF（对比图） ===
clear; close all; clc;

%% 1) 基本参数
Tp   = 2.0;        % 预定义时间
p    = 0.8;        % 时间归一化幂次（G-PPF 用）
a    = 0.20;       % 稳态下界
iota = 3.0;        % G-PPF 门控平滑
T_end = 4.0;       % 展示到 4 s

% G-PPF 呼吸/自适应参数（较强，确保扰动段包络）
Sigma_max  = 1.00;     % 最大放宽倍数：rho = a*(1 + g*Sigma)（总幅 ~ a*(1+Sigma_max)）
k_d        = 2.50;     % 扰动权重
tau_d_post = 0.15;     % 扰动低通
k_e        = 1.20;     % 误差权重
tau_e      = 0.07;     % 误差低通

% G-PPF 前段 sigma(t)（饱和→残差）
sigma0     = 1.00;  sigma_min = 0.60;  sigma_max = 1.60;
k_u        = 1.00;  tau_u_pre = 0.15;

%% 2) 时间轴
N  = 4800;  t0 = 1e-12;
t  = linspace(t0, T_end, N)';   dt = t(2)-t(1);
idx_pre  = t <  Tp;
idx_post = t >= Tp;

%% 3) 误差（初始更大 & 扰动为单段正弦）
E0 = 5.0*a;                           % 初始误差较大
e  = E0*exp(-1.1*t);
% 前段小振荡
e(idx_pre)  = e(idx_pre)  + 0.10*a*sin(22*t(idx_pre));
% 扰动窗（单段正弦）
t_on=2.40; t_off=3.30; alpha=10;
seg = 0.5*(tanh(alpha*(t - t_on)) - tanh(alpha*(t - t_off)));  % 平滑 0-1 窗
f_sin = 2.0;   A_sin = 1.05*a;
e(idx_post) = e(idx_post) + A_sin .* seg(idx_post) .* sin(2*pi*f_sin*(t(idx_post)-t_on));

%% 4) —— 普通性能函数 ρ_std(t)（用于对比）——
rho0   = 3*a;                         % 初始管径（可调小些以产生初值越界）
eps_Tp = 0.01;                        % t=Tp 时剩余比例 1%
lambda = -log(eps_Tp)/Tp;             % e^{-lambda*Tp} = eps_Tp
rho_std = (rho0 - a)*exp(-lambda*t) + a;  % 全程收敛到 a
% （注意：ρ_std(0)=rho0 为有限值，不具“全局初值可行性”）

%% 5) —— 全局自适应性能函数 ρ_g-PPF(t) —— 
% 时间归一化与核/窗
b   = (t./Tp).^p;
s   = 1 - 3*b.^2 + 2*b.^3;
phi = -log(b);
% 门控（G-PPF）：Tp 后平滑开启
g    = zeros(N,1);  g(idx_post) = 1 - exp( -iota*(t(idx_post)-Tp).^2 );

% 前段：饱和→残差→sigma(t)
u_max = 1.5;   u_cmd = 2.1 + 0.9*sin(14*t);
DeltaU_pre = max(0, abs(u_cmd) - u_max);
taper = double(idx_pre).*(1 - b.^2);
ru_pre  = zeros(N,1);
for k=2:N
    ru_pre(k) = ru_pre(k-1) + dt*( -ru_pre(k-1) + DeltaU_pre(k).*taper(k) )/tau_u_pre;
end
sigma_adapt = min(max(sigma0 + k_u*ru_pre, sigma_min), sigma_max);

% 后段：扰动 r_d 与误差 r_e
Distb_post = 1.6 * seg .* sin(2*pi*f_sin*(t - t_on)) .* idx_post; % 扰动（单段正弦）
rd_post = zeros(N,1);
for k=2:N
    rd_post(k) = rd_post(k-1) + dt*( -rd_post(k-1) + Distb_post(k) )/tau_d_post;
end
re_post = zeros(N,1);
for k=2:N
    re_post(k) = re_post(k-1) + dt*( -re_post(k-1) + abs(e(k))*idx_post(k) )/tau_e;
end
Sigma0 = min(max(k_d*abs(rd_post) + k_e*re_post, 0), Sigma_max);

% 扰动窗内直接全开：g_eff≈1
g_eff = max(g, seg);

% 初算 ρ_g-PPF
Sigma = Sigma0;
rho_gppf = a + sigma_adapt .* phi .* s;               % pre
rho_gppf(idx_post) = a.*(1 + g_eff(idx_post).*Sigma(idx_post));  % post

% 包络修正：确保扰动段 |e| <= ρ_g-PPF（不改变普通 PPF）
gap  = abs(e) - rho_gppf;
need = max(gap,0) ./ max(a*g_eff,1e-9);
add  = min(need, Sigma_max - Sigma);
Sigma = Sigma + add;
rho_gppf(idx_post) = a.*(1 + g_eff(idx_post).*Sigma(idx_post));

%% 6) 绘图（对比：普通 PPF vs 全局自适应 G-PPF）
col_std  = [0.00 0.45 0.74];  % 蓝：普通 PPF
col_gppf = [0.85 0.33 0.10];  % 橙：全局自适应 G-PPF
col_err  = [0.00 0.50 0.00];  % 绿：误差

figure('Color','w','Units','inches','Position',[1 1 7.2 3.2]);
ax = axes('Position',[0.12 0.20 0.84 0.72]); hold(ax,'on'); box(ax,'on');
set(ax,'FontName','Times New Roman','LineWidth',0.9,'TickDir','out',...
       'XMinorGrid','on','YMinorGrid','on','GridAlpha',0.16,'MinorGridAlpha',0.09);

% ±ρ_std（普通 PPF）
plot(t,  rho_std, '--','LineWidth',1.7,'Color',col_std,  'DisplayName','\pm\rho_{std}(t)（普通PPF）');
plot(t, -rho_std, '--','LineWidth',1.7,'Color',col_std,  'HandleVisibility','off');

% ±ρ_g-PPF（全局自适应）
plot(t,  rho_gppf, '-','LineWidth',2.1,'Color',col_gppf,'DisplayName','\pm\rho_{g\text{-}PPF}(t)（全局自适应）');
plot(t, -rho_gppf, '-','LineWidth',2.1,'Color',col_gppf,'HandleVisibility','off');

% e(t)
plot(t, e, '-.','LineWidth',1.8,'Color',col_err,'DisplayName','e(t)');

xlim([0, T_end]);   Ycap = 3.0;  ylim([-Ycap, Ycap]);

% T_p 参考线（不进图例）
hTp = xline(Tp,'k--','T_p','LabelVerticalAlignment','bottom','LabelOrientation','horizontal');
set(hTp,'HandleVisibility','off');

xlabel('时间 t (s)'); ylabel('幅值');
title('普通性能函数 vs 全局自适应性能函数（G-PPF）');
leg = legend('Location','northeast'); set(leg,'Box','off');

%% 7) 图内 ±∞ 箭头（针对 G-PPF 边界，避免与曲线重合）
x_inf = 0.10;  Y = ylim;
quiver(x_inf,  0.55*Y(2),  0,  0.35*Y(2), 0, 'Color',col_gppf,'LineWidth',1.2,'MaxHeadSize',0.9);
text(   x_inf+0.06, 0.90*Y(2), '+\infty','FontName','Times New Roman','FontSize',10,...
       'Color',col_gppf,'HorizontalAlignment','left','VerticalAlignment','bottom');
quiver(x_inf,  0.55*Y(1),  0,  0.35*Y(1), 0, 'Color',col_gppf,'LineWidth',1.2,'MaxHeadSize',0.9);
text(   x_inf+0.06, 0.90*Y(1), '-\infty','FontName','Times New Roman','FontSize',10,...
       'Color',col_gppf,'HorizontalAlignment','left','VerticalAlignment','top');
