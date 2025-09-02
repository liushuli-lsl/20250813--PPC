function [rho, rho_dot, aux] = gppf(t, e, DeltaU, d, cfg)
% GPPF: Global Prescribed-Performance Function (vector form, per-channel)
% Inputs:
%   t       : current time (scalar)
%   e       : n×1 (用于 re 通道的 |e|，建议传 kappa.*e)
%   DeltaU  : n×1 执行器饱和残差 Δu = u - sat(u)
%   d       : n×1 扰动幅值估计（例如 |d_hat|），没有可传 0
%   cfg     : 参数结构体（见控制器内 cfg1/cfg2）
%             必填: id, Tp, p, a, sigma0, sigma_min, sigma_max, iota, Sigma_max,
%                   k_u,k_d,k_e, use_lpf, tau_u,tau_d,tau_e
%             可选: feed_only=true  时，仅更新内部 LPF，不返回有效 rho
% Outputs:
%   rho, rho_dot : n×1
%   aux          : struct，含 ru,rd,re,g,Sigma,sigma（便于调试/记录）

% -------- persistent 状态：为每个 id 维护一组滤波器 --------
persistent S
if isempty(S)
    S = containers.Map('KeyType','double','ValueType','any');
end
id = cfg.id;

% 初始化该 id 的状态
if ~isKey(S,id)
    S(id) = struct('t_prev', t, ...
                   'ru', zeros(size(e)), ...
                   'rd', zeros(size(e)), ...
                   're', zeros(size(e)));
end
Si = S(id);

% 估计 dt（适配可变步长 ODE 求解器）
dt = max( t - Si.t_prev, 0 );

% ---- 驱动信号：低通或直接使用（由 cfg.use_lpf 决定）----
if isfield(cfg,'use_lpf') && cfg.use_lpf
    if dt > 0
        Si.ru = Si.ru + dt*( -Si.ru + abs(DeltaU) )/cfg.tau_u;
        Si.rd = Si.rd + dt*( -Si.rd + abs(d)      )/cfg.tau_d;
        Si.re = Si.re + dt*( -Si.re + abs(e)      )/cfg.tau_e;
    end
else
    Si.ru = abs(DeltaU);
    Sim.rd= abs(d);
    Si.re = abs(e);
end

% 如果仅 feed（回灌），直接更新状态并返回（不改变当步 rho）
if isfield(cfg,'feed_only') && cfg.feed_only
    Si.t_prev = t;
    S(id) = Si;
    rho = []; rho_dot = []; aux = [];
    return
end

% ---- 归一化时间核与窗 ----
b     = max( (t/cfg.Tp)^cfg.p, 1e-6 );       % 避免 log(0)
phi   = -log(b);                              % phi(b)
s     = 1 - 3*b^2 + 2*b^3;                    % s(b)
phip  = -1/b;                                 % dphi/db
sp    = -6*b + 6*b^2;                         % ds/db
b_dot = (cfg.p/cfg.Tp) * (t/cfg.Tp)^(cfg.p-1);

% ---- sigma (pre-phase) 与 g, Sigma (post-phase) ----
sigma = proj(cfg.sigma0 + cfg.k_u*Si.ru + cfg.k_d*Si.rd, cfg.sigma_min, cfg.sigma_max);

if t < cfg.Tp
    g = 0;
else
    g = 1 - exp( -cfg.iota*(t - cfg.Tp)^2 );
end

Sigma = proj(cfg.k_u*Si.ru + cfg.k_d*Si.rd + cfg.k_e*Si.re, 0, cfg.Sigma_max);

% ==== 计算 rho 与 rho_dot ====
n = numel(e);
rho     = zeros(n,1);
rho_dot = zeros(n,1);

if t < cfg.Tp
    % 预定义时间内： a + sigma*phi*s  （sigma 视为慢变，忽略 sigma_dot）
    rho     = cfg.a + sigma.*(phi*s);
    rho_dot =          sigma.*( phip*s + phi*sp )*b_dot;
else
    % 收敛后： a*(1 + g*Sigma)
    gp        = 2*cfg.iota*(t - cfg.Tp)*exp(-cfg.iota*(t - cfg.Tp)^2);
    % LPF 的解析时间导数（有 dt>0 时有效；否则近似为 0）
    if dt > 0 && isfield(cfg,'use_lpf') && cfg.use_lpf
        ru_dot = ( -Si.ru + abs(DeltaU) )/cfg.tau_u;
        rd_dot = ( -Si.rd + abs(d)      )/cfg.tau_d;
        re_dot = ( -Si.re + abs(e)      )/cfg.tau_e;
        Sigma_dot = cfg.k_u*ru_dot + cfg.k_d*rd_dot + cfg.k_e*re_dot;
    else
        Sigma_dot = zeros(n,1);
    end
    rho     = cfg.a*(1 + g*Sigma);
    rho_dot = cfg.a*( gp*Sigma + g.*Sigma_dot );
end

% 更新状态
Si.t_prev = t;
S(id) = Si;

% 辅助量
aux = struct('ru',Si.ru,'rd',Si.rd,'re',Si.re, ...
             'g',g,'Sigma',Sigma,'sigma',sigma);
end

% ----- helpers -----
function y = proj(x, xmin, xmax)
    y = min(max(x, xmin), xmax);
end
