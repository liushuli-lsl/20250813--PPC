function [rho, rho_dot, aux] = gppf(t, n, e, DeltaU, d, cfg)
% GPPF: Global Prescribed-Performance Function (vector form, per-channel)
% - 0<t<Tp:   rho = a + sigma * phi(b) * s(b)
% - t>=Tp:    rho = a * (1 + Sigma),  无门控灵敏缩放（回滞触发+积分记忆+投影）
% - 返回 aux.ratio = dot{rho}/rho（含数值保护）
%
% 必要 cfg 字段（若缺失自动赋默认值）：
%   id, Tp, p, a
%   sigma0, k_u, k_d, k_e, sigma_min, sigma_max, tau_u, tau_d, tau_e
%   % 下面为 t>=Tp 新增参数（均可用默认值）
%   u_d, theta_on, theta_off, eps_sw, ...
%   lambda_S, k_r, k_u_post, k_d_post, k_e_post, k_j_post, ...
%   k_I, alpha, Sigma_max, ...
%   d_th, d_scale, q_th, e_scale, j_th, j_scale, ...
%   gamma_b, gamma_q, k_burst
%
% 说明：
%   - s_bar 用 DeltaU 代理饱和强度: s̄ ≈ ReLU_eps(|DeltaU|)/u_d （无需 v）
%   - j_bar 默认禁用（=0），若你将 |v_dot| 提供到 cfg.vdot，可开启

% ------------------------- defaults -------------------------
cfg = gppf_defaults(cfg);

% ------------------ persistent per-id state ------------------
persistent S
if isempty(S)
    S = containers.Map('KeyType','double','ValueType','any');
end
id = cfg.id;
if ~isKey(S,id)
    S(id) = init_state(n, t);
end
Si = S(id);

% ---------------------- parse inputs -------------------------
e      = e(:);                 if numel(e)==1,      e = repmat(e, n,1);      end
if nargin<4 || isempty(DeltaU), DeltaU = zeros(n,1); do_feed_u = false; else
    DeltaU = DeltaU(:);        if numel(DeltaU)==1, DeltaU = repmat(DeltaU,n,1); end
    do_feed_u = true;
end
if nargin<5 || isempty(d),      d = zeros(n,1);     do_feed_d = false; else
    d = d(:);                  if numel(d)==1,      d = repmat(d, n,1);      end
    do_feed_d = true;
end

% --------------------- time / LPF update ---------------------
dt = max(t - Si.t_prev, 0);                 % robust dt for variable-step
% LPF for ru, rd, re
if dt>0
    if do_feed_u, Si.ru = Si.ru + dt*( -Si.ru + abs(DeltaU) )/cfg.tau_u; end
    if do_feed_d, Si.rd = Si.rd + dt*( -Si.rd + abs(d)      )/cfg.tau_d; end
                Si.re = Si.re + dt*( -Si.re + abs(e)        )/cfg.tau_e;
end

% 仅回灌模式
if isfield(cfg,'feed_only') && cfg.feed_only
    Si.t_prev = t; S(id)=Si; rho=[]; rho_dot=[]; aux=[]; return;
end

% ---------------------- pre-phase (t<Tp) ---------------------
t_eps = 1e-6;     b_min = 1e-6;
t_eff = max(t, t_eps);
b     = max( (t_eff/cfg.Tp)^cfg.p, b_min );
phi   = -log(b);
s     = 1 - 3*b^2 + 2*b^3;
phip  = -1/b;
sp    = -6*b + 6*b^2;
b_dot = (cfg.p/cfg.Tp) * (t_eff/cfg.Tp)^(cfg.p-1);

sigma = proj(cfg.sigma0 + cfg.k_u*Si.ru + cfg.k_d*Si.rd + cfg.k_e*Si.re, ...
             cfg.sigma_min, cfg.sigma_max);

rho     = zeros(n,1);
rho_dot = zeros(n,1);
ratio   = zeros(n,1);

if t < cfg.Tp
    % 0 < t < Tp
    rho     = cfg.a + sigma.*(phi*s);
    rho_dot =          sigma.*( phip*s + phi*sp )*b_dot;

    if t <= t_eps
        ratio = zeros(n,1);
        rho_dot = clamp(rho_dot, 1e6);
    else
        ratio = rho_dot ./ max(rho, 1e-12);
    end

else
% -------------------- post-phase (t>=Tp) ---------------------
% 切换瞬间施加 C1 初值：Sigma(Tp)=0, xi(Tp)=0
    if Si.just_switched==false && Si.phase==0
        Si.Sigma = zeros(n,1);
        Si.xi    = zeros(n,1);
        Si.prev_sbar = zeros(n,1);
        Si.just_switched = true;
        Si.phase = 1;  % 已进入后段
    end

    % --- 平滑ReLU + 回滞触发（用DeltaU代理饱和强度）---
    ReLU = @(x,eps) 0.5*(x + sqrt(x.^2 + eps.^2));
    H    = @(x,eps) 0.5*(1 + tanh(x./eps));

    s_bar = ReLU(abs(DeltaU) - 0, cfg.eps_sw) / max(cfg.u_d, 1e-12); % 以 |DeltaU| 近似饱和强度
    w_bar = ReLU(abs(d)      - cfg.d_th, cfg.eps_sw) / max(cfg.d_scale, 1e-12);
    q_bar = ReLU(abs(e)      - cfg.q_th, cfg.eps_sw) / max(cfg.e_scale, 1e-12);

    % 输入速率指标（可选）：若 cfg.vdot 存在，则使用；否则置 0
    if isfield(cfg,'vdot') && ~isempty(cfg.vdot)
        vdot = cfg.vdot(:);  if numel(vdot)==1, vdot = repmat(vdot,n,1); end
        j_bar = ReLU(abs(vdot) - cfg.j_th, cfg.eps_sw) / max(cfg.j_scale,1e-12);
    else
        j_bar = zeros(n,1);
    end

    chi  = H(s_bar - cfg.theta_on,  cfg.eps_sw);
    zeta = H(cfg.theta_off - s_bar, cfg.eps_sw);

    % --- 显著化：指数增益 + 阈越爆发 ---
    if dt>0
        sbar_dot = (s_bar - Si.prev_sbar)/dt;
    else
        sbar_dot = zeros(n,1);
    end
    Si.prev_sbar = s_bar;

    gain_boost = exp( cfg.gamma_b*s_bar + cfg.gamma_q*q_bar );
    burst      = cfg.k_burst * max(sbar_dot,0) .* chi;

    % --- 组合驱动 ---
    U = (cfg.k_u_post.*s_bar + cfg.k_d_post.*w_bar + cfg.k_e_post.*q_bar + cfg.k_j_post.*j_bar) ...
        .* gain_boost + burst;

    % --- ODE: Sigma, xi  （Euler；若你是 RK，可把 dSigma/dxi 带入求解器）
    if dt>0
        dSigma   = -cfg.lambda_S.*Si.Sigma - cfg.k_r.*zeta.*Si.Sigma + chi.*U + cfg.k_I.*Si.xi;
        dxi      = -cfg.alpha.*Si.xi + chi.*U;
        Si.Sigma = Si.Sigma + dt.*dSigma;
        Si.xi    = Si.xi    + dt.*dxi;
        Si.Sigma = proj(Si.Sigma, 0, cfg.Sigma_max);  % 投影
    else
        dSigma = zeros(n,1);
    end

    % --- rho 与导数 ---
    rho     = cfg.a*(1 + Si.Sigma);
    rho_dot = cfg.a * dSigma;                    % 近似 \dot{\rho}=a\dot{\Sigma}
    ratio   = rho_dot ./ max(rho, 1e-12);

    % 附：把 s_bar/w_bar/... 暴露出去便于调参与绘图
    aux_post = struct('s_bar',s_bar,'w_bar',w_bar,'q_bar',q_bar,'j_bar',j_bar, ...
                      'chi',chi,'zeta',zeta,'gain_boost',gain_boost,'burst',burst);
end

% ---------------------- state & outputs ----------------------
Si.t_prev = t;
if t < cfg.Tp, Si.phase=0; Si.just_switched=false; end
S(id) = Si;

aux = struct('ru',Si.ru,'rd',Si.rd,'re',Si.re, ...
             'Sigma',Si.Sigma,'xi',Si.xi,'sigma',sigma, ...
             'ratio',ratio);
if exist('aux_post','var')
    aux = catstruct(aux, aux_post);
end
end

% ========================= helpers ==========================
function S = init_state(n, t0)
S = struct();
S.t_prev = t0;
S.ru = zeros(n,1);
S.rd = zeros(n,1);
S.re = zeros(n,1);
S.Sigma = zeros(n,1);
S.xi = zeros(n,1);
S.prev_sbar = zeros(n,1);
S.just_switched = false;
S.phase = 0;  % 0: pre, 1: post
end

function y = proj(x, xmin, xmax)
y = min(max(x, xmin), xmax);
end

function y = clamp(x, cap)
y = max(min(x, cap), -cap);
end

function out = catstruct(a,b)
% concatenate two scalar structs
f = fieldnames(b);
out = a;
for i=1:numel(f), out.(f{i}) = b.(f{i}); end
end

function cfg = gppf_defaults(cfg)
% fill defaults if fields missing
req = {'id','Tp','p','a'};
for k=1:numel(req)
    if ~isfield(cfg,req{k}), error('cfg.%s is required', req{k}); end
end
def = struct( ...
 'sigma0',0.5,'k_u',0.8,'k_d',0.6,'k_e',0.4,...
 'sigma_min',0.3,'sigma_max',2.0,'tau_u',0.2,'tau_d',0.25,'tau_e',0.2, ...
 'u_d',1.0, 'theta_on',0.01,'theta_off',0.005,'eps_sw',1e-3, ...
 'lambda_S',4.0,'k_r',8.0,'k_u_post',6.0,'k_d_post',2.5,'k_e_post',1.5,'k_j_post',1.0, ...
 'k_I',2.0,'alpha',8.0,'Sigma_max',3.5, ...
 'd_th',0.12,'d_scale',0.30,'q_th',0.06,'e_scale',0.30,'j_th',0.8,'j_scale',4.0, ...
 'gamma_b',3.0,'gamma_q',1.5,'k_burst',1.2 ...
 );
f = fieldnames(def);
for i=1:numel(f)
    if ~isfield(cfg,f{i}), cfg.(f{i}) = def.(f{i}); end
end
end



% function [rho, rho_dot, aux] = gppf(t, n,e, DeltaU, d, cfg)
% % GPPF: Global Prescribed-Performance Function (vector form, per-channel)
% % 关键改动：t=0 数值保护 + aux.ratio = (dot{rho}/rho) 的稳健输出（t=0 置 0）
% 
% % -------- persistent 状态：为每个 id 维护一组滤波器 --------
% persistent S
% if isempty(S)
%     S = containers.Map('KeyType','double','ValueType','any');
% end
% id = cfg.id;
% 
% % 初始化该 id 的状态
% if ~isKey(S,id)
%     S(id) = struct('t_prev', t, ...
%                    'ru', zeros(size(e)), ...
%                    'rd', zeros(size(e)), ...
%                    're', zeros(size(e)));
% end
% Si = S(id);
% 
% % 估计 dt（适配可变步长 ODE 求解器）
% dt = max(t - Si.t_prev, 0);
% 
% % % ---- 驱动信号：低通或直接使用（由 cfg.use_lpf 决定）----
% % if isfield(cfg,'use_lpf') && cfg.use_lpf
% %     if dt > 0
% %         Si.ru = Si.ru + dt*( -Si.ru + abs(DeltaU) )/cfg.tau_u;
% %         Si.rd = Si.rd + dt*( -Si.rd + abs(d)      )/cfg.tau_d;
% %         Si.re = Si.re + dt*( -Si.re + abs(e)      )/cfg.tau_e;
% %     end
% % else
% %     Si.ru = abs(DeltaU);
% %     Si.rd = abs(d);          % 修正原来的 Sim.rd
% %     Si.re = abs(e);
% % end
% % ---- parse optional inputs (robust) ----
% % 1) DeltaU: [] => do not update ru this step
% do_feed_u = true;
% if nargin < 3 || isempty(DeltaU)
%     do_feed_u = false; 
%     DeltaU = zeros(n,1);         % 占位，保证尺寸；但不用于更新
% else
%     DeltaU = DeltaU(:);
%     if numel(DeltaU)==1, DeltaU = repmat(DeltaU,n,1); end
% end
% 
% % 2) d: [] => treat as zeros (可选是否也跳过更新)
% do_feed_d = true;
% if nargin < 4 || isempty(d)
%     do_feed_d = false;
%     d = zeros(n,1);
% else
%     d = d(:);
%     if numel(d)==1, d = repmat(d,n,1); end
% end
% 
% % ---- LPF updates (only when we actually feed) ----
% if do_feed_u
%     Si.ru = Si.ru + dt*( -Si.ru + abs(DeltaU) )/cfg.tau_u;
% end
% if do_feed_d
%     Si.rd = Si.rd + dt*( -Si.rd + abs(d) )/cfg.tau_d;
% end
% % 误差通道通常每步都可以更新（或同样做保护）
% Si.re = Si.re + dt*( -Si.re + abs(e) )/cfg.tau_e;
% 
% % 如果仅 feed（回灌），直接更新状态并返回（不改变当步 rho）
% if isfield(cfg,'feed_only') && cfg.feed_only
%     Si.t_prev = t;
%     S(id) = Si;
%     rho = []; rho_dot = []; aux = [];
%     return
% end
% 
% % ---- t=0 数值保护：用 t_eff 与 b_min ----
% t_eps = 1e-6;                 % 极小起跳时间
% b_min = 1e-6;                 % 归一化时间下界，避免 log(0)
% t_eff = max(t, t_eps);        % 用 t_eff 代替 t
% 
% % ---- 归一化时间核与窗 ----
% b     = max( (t_eff/cfg.Tp)^cfg.p, b_min);   % 避免 log(0)
% phi   = -log(b);                              % phi(b)
% s     = 1 - 3*b^2 + 2*b^3;                   % s(b)
% phip  = -1/b;                                 % dphi/db
% sp    = -6*b + 6*b^2;                         % ds/db
% b_dot = (cfg.p/cfg.Tp) * (t_eff/cfg.Tp)^(cfg.p-1);
% 
% % ---- sigma (pre-phase) 与 g, Sigma (post-phase) ----
% sigma = proj(cfg.sigma0 + cfg.k_u*Si.ru + cfg.k_d*Si.rd+cfg.k_e*Si.re, cfg.sigma_min, cfg.sigma_max);
% 
% if t < cfg.Tp
%     g = 0;
% else
%     g = 1 - exp( -cfg.iota*(t - cfg.Tp)^2 );
% end
% 
% Sigma = proj(cfg.k_u*Si.ru + cfg.k_d*Si.rd + cfg.k_e*Si.re, 0, cfg.Sigma_max);
% 
% % ==== 计算 rho 与 rho_dot ====
% n = numel(e);
% rho     = zeros(n,1);
% rho_dot = zeros(n,1);
% ratio   = zeros(n,1);  % aux.ratio = dot{rho}/rho（稳健实现）
% 
% if t < cfg.Tp
%     % 预定义时间内： a + sigma*phi*s  （sigma 视为慢变，忽略 sigma_dot）
%     rho     = cfg.a + sigma.*(phi*s);
%     rho_dot =          sigma.*( phip*s + phi*sp )*b_dot;
% 
%     % —— 在 t 非常小时的稳健处理：
%     if t <= t_eps
%         % dot{rho}/rho 的右极限为 0：直接置 0，避免把 1/t 型带入控制律
%         ratio = zeros(n,1);
%         % 可选：对 rho_dot 做幅值限制，防止日志或乘法放大（不影响稳定性）
%         rho_dot_cap = 1e6;   % 如需更严，可调小些
%         rho_dot = max(min(rho_dot, rho_dot_cap), -rho_dot_cap);
%     else
%         ratio = rho_dot ./ max(rho, 1e-12);
%     end
% else
%     % 收敛后： a*(1 + g*Sigma)
%     gp = 2*cfg.iota*(t - cfg.Tp) * exp(-cfg.iota*(t - cfg.Tp)^2);
% 
%     % LPF 的解析时间导数（有 dt>0 时有效；否则近似为 0）
%     if dt > 0 && isfield(cfg,'use_lpf') && cfg.use_lpf
%         ru_dot = ( -Si.ru + abs(DeltaU) )/cfg.tau_u;
%         rd_dot = ( -Si.rd + abs(d)      )/cfg.tau_d;
%         re_dot = ( -Si.re + abs(e)      )/cfg.tau_e;
%         Sigma_dot = cfg.k_u*ru_dot + cfg.k_d*rd_dot + cfg.k_e*re_dot
%     else
%         Sigma_dot = zeros(n,1);
%     end
% 
%     rho     = cfg.a*(1 + g*Sigma);
%     rho_dot = cfg.a*( gp*Sigma + g.*Sigma_dot );
%     ratio   = rho_dot ./ max(rho, 1e-12);
% end
% 
% % 更新状态
% Si.t_prev = t;
% S(id) = Si;
% 
% % 辅助量
% aux = struct('ru',Si.ru,'rd',Si.rd,'re',Si.re, ...
%              'g',g,'Sigma',Sigma,'sigma',sigma, ...
%              'ratio',ratio);   % 用于控制律的 \dot\rho/\rho 稳健补偿
% end
% 
% % ----- helpers -----
% function y = proj(x, xmin, xmax)
%     y = min(max(x, xmin), xmax);
% end
