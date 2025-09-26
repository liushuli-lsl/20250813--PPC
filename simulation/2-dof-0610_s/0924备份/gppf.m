% 
% 
% function [rho, drho_over_rho, aux] = gppf(t, n, e, DeltaU, cfg)
% % GPPF (concise) — Global Prescribed-Performance Function
% % 0<t<Tp : rho_i = a + sigma_i * phi(b) * s(b)
% % t>=Tp  : rho_i = a * (1 + Sigma_i),
% %          dSigma_i = -λ(Σ_i-Σreq_i) - k_r*off_i*Σ_i
% % Channel-wise scaling by |DeltaU_i| & |e_i|/rho_i.
% 
% %% ---- required & defaults ----
% assert(isfield(cfg,'id') && isfield(cfg,'Tp'),'cfg.id & cfg.Tp required');
% dflt = struct( ...
%   'p',0.3,'a',0.01, ...                 % shape & base radius
%   'sigma0',0.5,'k_u',2.5,'k_e',1.5, ... % pre: sigma gain
%   'sigma_min',0.2,'sigma_max',80.0,'tau_u',0.20,'tau_e',0.20, ...
%   'ks_tight',1.2,'ks_relax',1.8, ...    % tighten/relax via slack/press
%   'u_d',1.0,'theta_on',0.01,'theta_off',0.004, ... % post: deadzone (normed)
%   'lambda_S',4.0,'k_r',6.0,'Sigma_max',80, ...
%   'delta_margin',0.05,'gamma_u',1.2,'gamma_e',1.0, 'blend_lpf',0.25, ...
%   'rmax_sigma',5.0,'rmax_Sigma',20,'hold_min',0.12);
% cfg = fill_defaults(cfg,dflt);
% 
% e = e(:);
% if nargin<4 || isempty(DeltaU), DeltaU = zeros(n,1); end
% DeltaU = DeltaU(:); if numel(DeltaU)==1, DeltaU = repmat(DeltaU,n,1); end
% epsr = 1e-9;
% 
% %% ---- persistent per id ----
% persistent M
% if isempty(M), M = containers.Map('KeyType','double','ValueType','any'); end
% if ~isKey(M,cfg.id)
%     S = struct('t_prev',t,'ru',zeros(n,1),'re',zeros(n,1), ...
%                'Sigma',zeros(n,1),'hold',zeros(n,1), ...
%                'sigma_old',zeros(n,1),'Sigma_old',zeros(n,1), ...
%                'rho_prev', cfg.a + cfg.sigma0*5 );     % 用于占用率估计
%     M(cfg.id) = S;
% end
% S = M(cfg.id);
% 
% %% ---- time & LPFs ----
% dt = max(t - S.t_prev, 0);
% if dt>0
%     S.ru = S.ru + dt*( -S.ru + abs(DeltaU) )/cfg.tau_u; % 饱和残差低通
%     S.re = S.re + dt*( -S.re + abs(e)      )/cfg.tau_e; % 误差幅值低通
% end
% 
% %% ---- outputs ----
% rho = zeros(n,1); rho_dot = zeros(n,1); drho_over_rho = zeros(n,1);
% 
% %% ================= pre-phase: 0 < t < Tp =================
% t_eps = 1e-6; t_eff = max(t,t_eps);  b = max((t_eff/cfg.Tp)^cfg.p,1e-6);
% phi = -log(b); s = 1-3*b^2+2*b^3;  phip=-1/b; sp=-6*b+6*b^2;
% b_dot = (cfg.p/cfg.Tp)*(t_eff/cfg.Tp)^(cfg.p-1);
% 
% if t < cfg.Tp
%     % two-way sigma: slack tightens, press relaxes
%     eta   = abs(e)/max(cfg.a,epsr);                 % occupancy (pre段用a近似rho)
%     slack = max(0,1-eta);
%     press = max(0,eta-1+0.05);
%     sigma = proj(cfg.sigma0 + cfg.k_u*S.ru + cfg.k_e*S.re ...
%                  - cfg.ks_tight*slack + cfg.ks_relax*press, ...
%                  cfg.sigma_min, cfg.sigma_max);
%     sigma = rate_limit_vec(sigma, S.sigma_old, dt, cfg.rmax_sigma);
%     S.sigma_old = sigma;
% 
%     rho     = cfg.a + sigma.*(phi*s);
%     rho_dot =          sigma.*(phip*s + phi*sp)*b_dot;
%     if t<=t_eps, drho_over_rho = zeros(n,1); rho_dot = clamp(rho_dot,1e6);
%     else,        drho_over_rho = rho_dot ./ max(rho,1e-12);
%     end
% 
% else
% %% ================= post-phase: t >= Tp =================
% % --- 1) 误差占用率（用上一步 rho 以避免代数环）---
% rho_guard = max(S.rho_prev, cfg.a);              % n×1
% occ = abs(e) ./ max(rho_guard, epsr);            % |e|/rho  (占用率)
% % 安全裕度阈：occ 超过 (1 - delta) 触发放宽；加入 on/off 双阈值避免抖振
% m_on  = occ - (1 - cfg.delta_margin) - cfg.theta_on;
% m_off = occ - (1 - cfg.delta_margin) - cfg.theta_off;
% m_on  = max(0, m_on);
% m_off = max(0, m_off);
% 
% % --- 2) 饱和强度（归一化 + on/off）---
% u_on  = S.ru/cfg.u_d - cfg.theta_on;
% u_off = S.ru/cfg.u_d - cfg.theta_off;
% u_on  = max(0, u_on);
% u_off = max(0, u_off);
% 
% % --- 3) 组合得到需求 Σ (ReLU 型)，并裁剪到 [0, Sigma_max] ---
% Sigma_req = cfg.gamma_u * u_on + cfg.gamma_e * m_on;        % 放宽需求
% Sigma_req = proj(Sigma_req, 0, cfg.Sigma_max);
% 
% % --- 4) 计算 off 指示：两者都低于 off 阈时，增加回收项 ---
% off_flag = (u_off<=0) & (m_off<=0);               % 逻辑矩阵
% off = double(off_flag);
% 
% % --- 5) 连续时间更新：dΣ = -λ(Σ-Σreq) - k_r*off*Σ  ---
% if dt > 0
%     dSigma = -cfg.lambda_S*(S.Sigma - Sigma_req) - cfg.k_r.*off.*S.Sigma;
%     % 速率限制，防止数值抖动
%     dSigma = min(max(dSigma, -cfg.rmax_Sigma), cfg.rmax_Sigma);
%     S.Sigma = proj(S.Sigma + dSigma*dt, 0, cfg.Sigma_max);
% else
%     dSigma = zeros(n,1);
% end
% 
% % --- 6) 输出与导数 ---
% a = cfg.a;
% rho           = a*(1+ S.Sigma);
% rho_dot       = a * dSigma;                      % 连续时间表达
% drho_over_rho = rho_dot ./ max(rho, 1e-12);
% 
% % 诊断输出（可选）
% aux_post = struct('occ',occ,'Sigma_req',Sigma_req,'off',off,'ru',S.ru,'re',S.re);
% end
% 
% %% ---- persist & aux ----
% S.t_prev   = t;
% S.rho_prev = rho;          % 存储本步 rho 供下步占用率估计
% M(cfg.id)  = S;
% 
% aux = struct('Sigma',S.Sigma);
% if exist('Sigma_req','var'), aux.Sigma_req = Sigma_req; end
% if exist('aux_post','var'),  aux.post = aux_post; end
% 
% end
% % ================= helpers =================
% function y = proj(x,xmin,xmax), y = min(max(x,xmin),xmax); end
% function y = clamp(x,cap),      y = max(min(x,cap),-cap);  end
% function y = rate_limit_vec(xn,xo,dt,rmax)
% if isempty(xo) || dt<=0, y = xn; else
%     dx = xn - xo; r = rmax*dt; dx = min(max(dx,-r),r); y = xo + dx;
% end
% end
% function cfg = fill_defaults(cfg, d)
% f = fieldnames(d);
% for k=1:numel(f), if ~isfield(cfg,f{k}), cfg.(f{k}) = d.(f{k}); end, end
% end
% 
% 
% 










function [rho, rho_dot, aux] = gppf(t, n,e, DeltaU, d, cfg)
% GPPF: Global Prescribed-Performance Function (vector form, per-channel)
% 关键改动：t=0 数值保护 + aux.ratio = (dot{rho}/rho) 的稳健输出（t=0 置 0）
% d=zeros(n,1);
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
dt = max(t - Si.t_prev, 0);

do_feed_u = true;
if nargin < 3 || isempty(DeltaU)
    do_feed_u = false; 
    DeltaU = zeros(n,1);         % 占位，保证尺寸；但不用于更新
else
    DeltaU = DeltaU(:);
    if numel(DeltaU)==1, DeltaU = repmat(DeltaU,n,1); end
end

% 2) d: [] => treat as zeros (可选是否也跳过更新)
do_feed_d = true;
if nargin < 4 || isempty(d)
    do_feed_d = false;
    d = zeros(n,1);
else
    d = d(:);
    if numel(d)==1, d = repmat(d,n,1); end
end
DeltaU
d
% ---- LPF updates (only when we actually feed) ----
if do_feed_u
    Si.ru = Si.ru + dt*( -Si.ru + abs(DeltaU) )/cfg.tau_u;
end
if do_feed_d
    Si.rd = Si.rd + dt*( -Si.rd + abs(d) )/cfg.tau_d;
end
% 误差通道通常每步都可以更新（或同样做保护）
Si.re = Si.re + dt*( -Si.re + abs(e) )/cfg.tau_e;

% 如果仅 feed（回灌），直接更新状态并返回（不改变当步 rho）
if isfield(cfg,'feed_only') && cfg.feed_only
    Si.t_prev = t;
    S(id) = Si;
    rho = []; rho_dot = []; aux = [];
    return
end

% ---- t=0 数值保护：用 t_eff 与 b_min ----
t_eps = 1e-6;                 % 极小起跳时间
b_min = 1e-6;                 % 归一化时间下界，避免 log(0)
t_eff = max(t, t_eps);        % 用 t_eff 代替 t

% ---- 归一化时间核与窗 ----
b     = max( (t_eff/cfg.Tp)^cfg.p, b_min);   % 避免 log(0)
phi   = -log(b);                              % phi(b)
s     = 1 - 3*b^2 + 2*b^3;                   % s(b)
phip  = -1/b;                                 % dphi/db
sp    = -6*b + 6*b^2;                         % ds/db
b_dot = (cfg.p/cfg.Tp) * (t_eff/cfg.Tp)^(cfg.p-1);

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
ratio   = zeros(n,1);  % aux.ratio = dot{rho}/rho（稳健实现）

if t < cfg.Tp
    % 预定义时间内： a + sigma*phi*s  （sigma 视为慢变，忽略 sigma_dot）
    rho     = cfg.a + sigma.*(phi*s);
    rho_dot =          sigma.*( phip*s + phi*sp )*b_dot;

    % —— 在 t 非常小时的稳健处理：
    if t <= t_eps
        % dot{rho}/rho 的右极限为 0：直接置 0，避免把 1/t 型带入控制律
        ratio = zeros(n,1);
        % 可选：对 rho_dot 做幅值限制，防止日志或乘法放大（不影响稳定性）
        rho_dot_cap = 1e6;   % 如需更严，可调小些
        rho_dot = max(min(rho_dot, rho_dot_cap), -rho_dot_cap);
    else
        ratio = rho_dot ./ max(rho, 1e-12);
    end
else
    % 收敛后： a*(1 + g*Sigma)
    gp = 2*cfg.iota*(t - cfg.Tp) * exp(-cfg.iota*(t - cfg.Tp)^2);

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
    ratio   = rho_dot ./ max(rho, 1e-12);
end

% 更新状态
Si.t_prev = t;
S(id) = Si;

% 辅助量
aux = struct('ru',Si.ru,'rd',Si.rd,'re',Si.re, ...
             'g',g,'Sigma',Sigma,'sigma',sigma, ...
             'ratio',ratio);   % 用于控制律的 \dot\rho/\rho 稳健补偿
end

% ----- helpers -----
function y = proj(x, xmin, xmax)
    y = min(max(x, xmin), xmax);
end










% 
% function [rho, drho_over_rho, aux] = gppf(t, n, e, DeltaU, d, cfg)
% % GPPF  Global Prescribed-Performance Function （兼容调用：[rho, drho] = gppf(t,n,e,DeltaU,d,cfg)）
% % 输出：
% %   rho            : n×1 性能半径
% %   drho_over_rho  : n×1 比值 \dot\rho / \rho   （用于 BLF 抵消项）
% %   aux            : 诊断结构体（可选）
% %
% % 说明：
% % - 0<t<Tp:   rho = a + sigma * phi(b) * s(b)
% % - t>=Tp :   rho = a * (1 + Sigma)  （无门控灵敏缩放：回滞触发+积分+投影）
% % - 支持多实例：cfg.id 为实例 id（double）
% %
% % 建议：第二个输出就是你代码里的 drho1
% 
% %% --------------- 默认参数检查与补全 ---------------
% if ~isfield(cfg,'id');       error('cfg.id 必填'); end
% if ~isfield(cfg,'Tp');       error('cfg.Tp 必填'); end
% if ~isfield(cfg,'p');        cfg.p = 1.0; end
% if ~isfield(cfg,'a');        cfg.a = 0.12; end
% 
% d=zeros(n,1)
% 
% % pre-Tp
% def_sigma = struct('sigma0',0.5,'k_u',2.5,'k_d',0.6,'k_e',2.5, ...
%     'sigma_min',0.3,'sigma_max',100,'tau_u',0.2,'tau_d',0.25,'tau_e',0.2);
% cfg = fill_defaults(cfg, def_sigma);
% 
% % post-Tp（灵敏缩放）
% def_post = struct( ...
%     'u_d',1.0,'theta_on',0.01,'theta_off',0.005,'eps_sw',1e-3, ...
%     'lambda_S',4.0,'k_r',8.0, ...
%     'k_u_post',6.0,'k_d_post',2.5,'k_e_post',1.5,'k_j_post',1.0, ...
%     'k_I',2.0,'alpha',8.0,'Sigma_max',3.5, ...
%     'd_th',0.12,'d_scale',0.30,'q_th',0.06,'e_scale',0.30,'j_th',0.8,'j_scale',4.0, ...
%     'gamma_b',3.0,'gamma_q',1.5,'k_burst',0.8);
% cfg = fill_defaults(cfg, def_post);
% 
% %% --------------- 形参规整 ---------------
% e = e(:);
% if nargin < 4 || isempty(DeltaU), DeltaU = zeros(n,1); do_feed_u = false; else
%     DeltaU = DeltaU(:); if numel(DeltaU)==1, DeltaU = repmat(DeltaU,n,1); end
%     do_feed_u = true;
% end
% if nargin < 5 || isempty(d), d = zeros(n,1); do_feed_d = false; else
%     d = d(:); if numel(d)==1, d = repmat(d,n,1); end
%     do_feed_d = true;
% end
% 
% %% --------------- 持久化状态（每个 id 一份） ---------------
% persistent MAP
% if isempty(MAP)
%     MAP = containers.Map('KeyType','double','ValueType','any');
% end
% id = cfg.id;
% if ~isKey(MAP,id)
%     MAP(id) = init_state(n, t);
% end
% S = MAP(id);
% 
% %% --------------- 时间步与 LPF 更新 ---------------
% dt = max(t - S.t_prev, 0);
% if dt > 0
%     if do_feed_u
%         S.ru = S.ru + dt*( -S.ru + abs(DeltaU) )/cfg.tau_u;
%     end
%     if do_feed_d
%         S.rd = S.rd + dt*( -S.rd + abs(d) )/cfg.tau_d;
%     end
%     S.re = S.re + dt*( -S.re + abs(e) )/cfg.tau_e;
% end
% 
% %% --------------- pre-phase: 0 < t < Tp ---------------
% t_eps = 1e-6; b_min = 1e-6;
% t_eff = max(t, t_eps);
% 
% b    = max( (t_eff/cfg.Tp)^cfg.p, b_min );
% phi  = -log(b);
% s    = 1 - 3*b^2 + 2*b^3;
% phip = -1/b;
% sp   = -6*b + 6*b^2;
% b_dot= (cfg.p/cfg.Tp) * (t_eff/cfg.Tp)^(cfg.p-1);
% 
% sigma = proj(cfg.sigma0 + cfg.k_u*S.ru + cfg.k_d*S.rd + cfg.k_e*S.re, ...
%              cfg.sigma_min, cfg.sigma_max);
% 
% rho           = zeros(n,1);
% rho_dot       = zeros(n,1);
% drho_over_rho = zeros(n,1);
% 
% if t < cfg.Tp
%     rho     = cfg.a + sigma.*(phi*s);
%     rho_dot =          sigma.*( phip*s + phi*sp )*b_dot;
% 
%     if t <= t_eps
%         drho_over_rho = zeros(n,1);
%         rho_dot = clamp(rho_dot, 1e6);
%     else
%         drho_over_rho = rho_dot ./ max(rho, 1e-12);
%     end
% 
% else
% %% --------------- post-phase: t >= Tp  （无门控灵敏缩放） ---------------
%  % ====== 指标（更稳健地“看见”饱和） ======
% ReLU = @(x,eps) 0.5*(x + sqrt(x.^2 + eps.^2));
% H    = @(x,eps) 0.5*(1 + tanh(x./eps));
% 
% s_from_DU = ReLU(abs(DeltaU), cfg.eps_sw) / max(cfg.u_d,1e-12);
% s_from_v  = zeros(n,1);
% if isfield(cfg,'vabs') && ~isempty(cfg.vabs)
%     s_from_v = ReLU(max(cfg.vabs - cfg.u_d, 0), cfg.eps_sw) / max(cfg.u_d,1e-12);
% end
% s_from_u  = zeros(n,1);
% if isfield(cfg,'uabs') && ~isempty(cfg.uabs)
%     s_from_u = ReLU(max(cfg.uabs - cfg.u_d, 0), cfg.eps_sw) / max(cfg.u_d,1e-12);
% end
% s_bar = max(s_from_DU, max(s_from_v, s_from_u));
% 
% % —— 误差自适应尺度 —— %
% if ~isfield(S,'e_rms'), S.e_rms = 1e-3*ones(n,1); end
% tau_es = 0.15;
% if dt>0, S.e_rms = S.e_rms + dt*( -S.e_rms + abs(e) )/tau_es; end
% e_scale_adapt = max(2.0*S.e_rms, 1e-4);
% q_bar = ReLU(abs(e), cfg.eps_sw) ./ e_scale_adapt;
% 
% % —— 触发/释放（更灵敏+最小占空） —— %
% if ~isfield(cfg,'theta_on'),  cfg.theta_on  = 0.005; end
% if ~isfield(cfg,'theta_off'), cfg.theta_off = 0.002; end
% chi_u = H(s_bar - cfg.theta_on,      cfg.eps_sw);
% chi_q = H(q_bar - 0.5*cfg.theta_on,  cfg.eps_sw);
% chi   = max(chi_u, chi_q);
% chi   = max(chi, 0.10);
% zeta  = H(cfg.theta_off - s_bar, cfg.eps_sw);
% 
% % —— 显著化 & 阈越爆发 —— %
% if dt>0, sbar_dot = (s_bar - S.prev_sbar)/dt; else, sbar_dot = zeros(n,1); end
% S.prev_sbar = s_bar;
% if ~isfield(cfg,'gamma_b'), cfg.gamma_b = 2.0; end
% if ~isfield(cfg,'gamma_q'), cfg.gamma_q = 1.0; end
% if ~isfield(cfg,'k_burst'), cfg.k_burst = 0.8; end
% gain_boost  = exp( cfg.gamma_b*s_bar + cfg.gamma_q*q_bar );
% burst       = cfg.k_burst * max(sbar_dot,0) .* chi;
% 
% % —— 组合驱动（仅饱和+误差；加基线） —— %
% if ~isfield(cfg,'k_u_post'), cfg.k_u_post = 10.0; end
% if ~isfield(cfg,'k_e_post'), cfg.k_e_post =  2.0; end
% U_base = 0.10 * q_bar;
% U = (cfg.k_u_post.*s_bar + cfg.k_e_post.*q_bar) .* gain_boost + burst + U_base;
% 
% % ODE、rho、drho_over_rho 保持你原先写法（只用 U、chi、zeta 更新）
% 
% 
% % ====== ODE: Σ, ξ ======
% if dt > 0
%     dSigma   = -cfg.lambda_S.*S.Sigma - cfg.k_r.*zeta.*S.Sigma + chi.*U + cfg.k_I.*S.xi;
%     dxi      = -cfg.alpha.*S.xi + chi.*U;
%     % 稍放宽上限，避免太早顶死（可按需调整）
%     S.Sigma  = proj(S.Sigma + dt.*dSigma, 0, max(cfg.Sigma_max, 3.0));
%     S.xi     = S.xi + dt.*dxi;
% else
%     dSigma = zeros(n,1);
% end
% 
% % ====== rho 与导数比 ======
% rho           = cfg.a*(1 + S.Sigma);
% rho_dot       = cfg.a * dSigma;                        % 近似 a·Σ̇
% drho_over_rho = rho_dot ./ max(rho, 1e-12);
% 
% % ====== 便于调参的辅助量 ======
% aux_post.s_bar = s_bar;
% aux_post.q_bar = q_bar;
% aux_post.chi   = chi;
% aux_post.zeta  = zeta;
% aux_post.gain_boost = gain_boost;
% aux_post.burst      = burst;
% 
% 
% if t > cfg.Tp && (mod(round(t*1000), 200)==0)  % 每0.2s打印一次
%     fprintf('t=%.2f | maxΔU=%.3g s_bar=%.3g | max|e|=%.3g q_bar=%.3g | chi=%.2f | Sigma=%.2f\n', ...
%         t, max(abs(DeltaU)), max(abs(e)));
% end
% 
% end
% 
% %% --------------- 更新状态并输出 ---------------
% S.t_prev = t;
% if t < cfg.Tp
%     S.phase = 0;
% end
% MAP(id) = S;
% 
% aux = struct();
% aux.ru = S.ru; aux.rd = S.rd; aux.re = S.re;
% aux.sigma = sigma; aux.Sigma = S.Sigma; aux.xi = S.xi;
% aux.ratio = drho_over_rho;
% if exist('aux_post','var')
%     % 手动合并，避免 catstruct 依赖
%     fields = fieldnames(aux_post);
%     for i = 1:numel(fields)
%         aux.(fields{i}) = aux_post.(fields{i});
%     end
% end
% 
% end % ===== END of gppf =====
% 
% %% -------------------- 本地小工具 --------------------
% function y = proj(x, xmin, xmax)
% y = min(max(x, xmin), xmax);
% end
% 
% function y = clamp(x, cap)
% y = max(min(x, cap), -cap);
% end
% 
% function S = init_state(n, t0)
% S = struct();
% S.t_prev = t0;
% S.ru = zeros(n,1);
% S.rd = zeros(n,1);
% S.re = zeros(n,1);
% S.Sigma = zeros(n,1);
% S.xi = zeros(n,1);
% S.prev_sbar = zeros(n,1);
% S.phase = 0; % 0: pre, 1: post
% end
% 
% function cfg = fill_defaults(cfg, def)
% f = fieldnames(def);
% for i = 1:numel(f)
%     k = f{i};
%     if ~isfield(cfg,k)
%         cfg.(k) = def.(k);
%     end
% end
% end



