% function [rho_p,rho_m,rho, rho_dot, aux] = gppf(t, n, e, DeltaU, d, cfg)
% % GPPF with Asymmetric/Tunnel Envelopes (Scheme A)
% % - Backward compatible: outputs rho, rho_dot (symmetric) unchanged.
% % - New (in aux.tunnel): rho_p/rho_m, rho_dot_p/rho_dot_m, ratio_p/ratio_m,
% %   delta, beta, s, zeta for tunnel (asymmetric) envelopes.
% %
% % Required cfg fields (existing + new):
% %   Tp, p, a, iota, sigma0, sigma_min, sigma_max, k_u, k_d, k_e, ...
% %   tau_u, tau_d, tau_e, Sigma_max,
% %   % NEW for tunnel:
% %   delta0, delta_max, c_u, c_d, c_e, beta0, gamma_beta
% %
% % Example defaults (建议):
% %   cfg.delta0=0.1; cfg.delta_max=0.5; cfg.c_u=1; cfg.c_d=1; cfg.c_e=1;
% %   cfg.beta0=0.3*cfg.a; cfg.gamma_beta=3;
% % —— 新增（方案A：隧道/非对称）——
% cfg.delta0     = 0.10;        % 基础相对涨宽（0~delta_max）
% cfg.delta_max  = 0.50;        % 相对涨宽上限（建议 0.3~0.6）
% cfg.c_u        = 1.0;         % ru -> delta 的权重（与 k_u 同量级）
% cfg.c_d        = 0.5;         % rd -> delta 的权重
% cfg.c_e        = 1.0;         % re -> delta 的权重
% cfg.beta0      = 0.30*cfg.a;  % 方向偏置幅度（与 a 同量级）
% cfg.gamma_beta = 3.0;         % 偏置方向灵敏度（tanh 的斜率）
% % pre-phase (t < T_p) 调整旋钮
% cfg.pre_q           = 2.0;      % w_pre(b) = (1-b)^q，q越大，越偏向早期更放宽
% cfg.sigma_max_pre   = 3.0;      % 仅在 t<Tp 生效的更大上限（> sigma_max）
% cfg.pre_gain        = 1.0;      % 对 ru/rd/re 的额外灵敏度（只在 t<Tp）
% cfg.pre_boost_gain  = 0.6;      % “脉冲放宽”幅度（只在 t<Tp）
% % -------- persistent 状态：为每个 id 维护一组滤波器 --------
% persistent S
% if isempty(S), S = containers.Map('KeyType','double','ValueType','any'); end
% id = cfg.id;
% 
% % 初始化该 id 的状态
% if ~isKey(S,id)
%     S(id) = struct('t_prev', t, ...
%                    'ru', zeros(n,1), ...
%                    'rd', zeros(n,1), ...
%                    're', zeros(n,1));
% end
% Si = S(id);
% 
% % ------- 输入与维度保护 -------
% if nargin < 3 || isempty(e), e = zeros(n,1); else, e = e(:); end
% 
% do_feed_u = true;
% if nargin < 4 || isempty(DeltaU)
%     do_feed_u = false;
%     DeltaU = zeros(n,1);
% else
%     DeltaU = DeltaU(:);
%     if numel(DeltaU)==1, DeltaU = repmat(DeltaU,n,1); end
% end
% 
% do_feed_d = true;
% if nargin < 5 || isempty(d)
%     do_feed_d = false;
%     d = zeros(n,1);
% else
%     d = d(:);
%     if numel(d)==1, d = repmat(d,n,1); end
% end
% 
% % ------- 估计 dt（适配变步长 ODE） -------
% dt = max(t - Si.t_prev, 0);
% 
% % ------- 一阶LPF状态更新（仅在feed时更新） -------
% if do_feed_u
%     Si.ru = Si.ru + dt*( -Si.ru + abs(DeltaU) )/cfg.tau_u;
% end
% if do_feed_d
%     Si.rd = Si.rd + dt*( -Si.rd + abs(d) )/cfg.tau_d;
% end
% % 误差通道通常每步更新
% Si.re = Si.re + dt*( -Si.re + abs(e) )/cfg.tau_e;
% 
% % 若仅回灌
% if isfield(cfg,'feed_only') && cfg.feed_only
%     Si.t_prev = t; S(id) = Si;
%     rho = []; rho_dot = []; aux = [];
%     return
% end
% 
% % ------- t=0 数值保护 -------
% t_eps = 1e-6;                % 极小起跳时间
% b_min = 1e-6;                % 避免 log(0)
% t_eff = max(t, t_eps);
% 
% % ------- 归一化时间核与窗 -------
% b     = max( (t_eff/cfg.Tp)^cfg.p, b_min);   % b in (b_min,1]
% phi   = -log(b);                              % phi(b)
% s     = 1 - 3*b^2 + 2*b^3;                    % s(b)
% phip  = -1/b;                                 % dphi/db
% sp    = -6*b + 6*b^2;                         % ds/db
% b_dot = (cfg.p/cfg.Tp) * (t_eff/cfg.Tp)^(cfg.p-1);
% 
% % ------- sigma, g, Sigma -------
% sigma = proj(cfg.sigma0 + cfg.k_u*Si.ru + cfg.k_d*Si.rd + cfg.k_e*Si.re, ...
%              cfg.sigma_min, cfg.sigma_max);
% 
% if t < cfg.Tp
%     g = 0;
% else
%     g = 1 - exp( -cfg.iota*(t - cfg.Tp)^2 );  % C^1 门控
% end
% 
% Sigma = proj(cfg.k_u*Si.ru + cfg.k_d*Si.rd + cfg.k_e*Si.re, 0, cfg.Sigma_max);
% 
% % ------- 计算对称 rho 与 rho_dot（向后兼容） -------
% rho     = zeros(n,1);
% rho_dot = zeros(n,1);
% ratio   = zeros(n,1);  % dot{rho}/rho（稳健实现）
% 
% if t < cfg.Tp
%     rho     = cfg.a + sigma.*(phi*s);
%     rho_dot =          sigma.*( phip*s + phi*sp )*b_dot;
% 
%     if t <= t_eps
%         ratio = zeros(n,1);
%         rho_dot_cap = 1e6;
%         rho_dot = max(min(rho_dot, rho_dot_cap), -rho_dot_cap);
%     else
%         ratio = rho_dot ./ max(rho, 1e-12);
%     end
% else
%     gp = 2*cfg.iota*(t - cfg.Tp) * exp(-cfg.iota*(t - cfg.Tp)^2);
% 
%     if dt > 0
%         ru_dot = ( -Si.ru + abs(DeltaU) )/cfg.tau_u;
%         rd_dot = ( -Si.rd + abs(d)      )/cfg.tau_d;
%         re_dot = ( -Si.re + abs(e)      )/cfg.tau_e;
%         Sigma_dot_raw = cfg.k_u*ru_dot + cfg.k_d*rd_dot + cfg.k_e*re_dot;
%         % Sigma_dot 经过投影梯度
%         Sigma_dot = proj_grad(Sigma_dot_raw, ...
%                               cfg.k_u*ru_dot + cfg.k_d*rd_dot + cfg.k_e*re_dot, ...
%                               0, cfg.Sigma_max, Sigma);
%     else
%         Sigma_dot = zeros(n,1);
%     end
% 
%     rho     = cfg.a*(1 + g*Sigma);
%     rho_dot = cfg.a*( gp*Sigma + g.*Sigma_dot );
%     ratio   = rho_dot ./ max(rho, 1e-12);
% end
% 
% % ================== 方案A：非对称/隧道 ==================
% % delta(t) （相对涨宽）与 beta(t)（方向偏置）
% delta_raw = cfg.delta0 + cfg.c_u*Si.ru + cfg.c_d*Si.rd + cfg.c_e*Si.re;
% delta     = proj(delta_raw, 0, cfg.delta_max);
% 
% % beta(t) = beta0 * g(t) * tanh(gamma_beta * e)
% beta      = cfg.beta0 .* g .* tanh(cfg.gamma_beta .* e);
% 
% % 导数：delta_dot, beta_dot（用于 BLF 中的 \dot\rho^\pm / \rho^\pm）
% if dt > 0
%     ru_dot = ( -Si.ru + abs(DeltaU) )/cfg.tau_u;
%     rd_dot = ( -Si.rd + abs(d)      )/cfg.tau_d;
%     re_dot = ( -Si.re + abs(e)      )/cfg.tau_e;
% 
%     delta_dot_raw = cfg.c_u*ru_dot + cfg.c_d*rd_dot + cfg.c_e*re_dot;
%     delta_dot = proj_grad(delta_dot_raw, delta_dot_raw, 0, cfg.delta_max, delta);
% 
%     gp = 2*cfg.iota*(t - cfg.Tp) * exp(-cfg.iota*(t - cfg.Tp)^2); % g'(t)
%     beta_dot = cfg.beta0 .* ( gp .* tanh(cfg.gamma_beta.*e) ...
%                + g .* (cfg.gamma_beta .* sech2(cfg.gamma_beta.*e)) .* e_dot_placeholder(n) );
%     % 说明：e_dot_placeholder 为站位项。若你在此函数中没有 e 的动态，
%     %       可置零（即只保留门控项 gp*tanh(...)），或在调用处传入 \dot e。
%     beta_dot = beta_dot; %#ok<NASGU>
% else
%     delta_dot = zeros(n,1);
%     beta_dot  = zeros(n,1);
% end
% 
% % 上/下隧道边界及其导数
% rho_p     = rho .* (1 + delta) + beta;
% rho_m     = rho .* (1 - delta) - beta;
% 
% % \dot\rho^\pm = \dot\rho(1\pm\delta) \pm rho*delta_dot \pm beta_dot
% rho_dot_p = rho_dot .* (1 + delta) + rho .* delta_dot + beta_dot;
% rho_dot_m = rho_dot .* (1 - delta) - rho .* delta_dot - beta_dot;
% 
% % 稳健比率（避免 0 除）
% ratio_p = rho_dot_p ./ max(rho_p, 1e-12);
% ratio_m = rho_dot_m ./ max(rho_m, 1e-12);
% 
% % 隧道归一化映射 s \in (0,1) 与 zeta
% % s = (e + rho_m) / (rho_p + rho_m)
% den_tun = max(rho_p + rho_m, 1e-12);
% s_norm  = (e + rho_m) ./ den_tun;
% s_norm  = min(max(s_norm, 1e-12), 1-1e-12);   % 数值保护
% zeta    = log( s_norm ./ (1 - s_norm) );
% 
% % ================== 保存与输出 ==================
% Si.t_prev = t;
% S(id) = Si;
% 
% aux = struct();
% aux.ru = Si.ru; aux.rd = Si.rd; aux.re = Si.re;
% aux.g = g; aux.Sigma = Sigma; aux.sigma = sigma; aux.ratio = ratio;
% 
% % 隧道相关输出
% aux.tunnel = struct();
% aux.tunnel.delta = delta;
% aux.tunnel.beta  = beta;
% aux.tunnel.rho_p = rho_p;
% aux.tunnel.rho_m = rho_m;
% aux.tunnel.rho_dot_p = rho_dot_p;
% aux.tunnel.rho_dot_m = rho_dot_m;
% aux.tunnel.ratio_p = ratio_p;
% aux.tunnel.ratio_m = ratio_m;
% aux.tunnel.s = s_norm;
% aux.tunnel.zeta = zeta;
% 
% end
% 
% % ---------- helpers ----------
% function y = proj(x, xmin, xmax)
%     y = min(max(x, xmin), xmax);
% end
% 
% % 投影梯度：当变量在区间内部时返回 raw；在边界且梯度“指向外”时置零；否则保留。
% function y = proj_grad(raw, grad, xmin, xmax, x)
%     y = grad;
%     mask_low  = (x <= xmin + 1e-12) & (grad < 0);
%     mask_high = (x >= xmax - 1e-12) & (grad > 0);
%     y(mask_low | mask_high) = 0;
% end
% 
% function y = sech2(x)
%     % sech^2(x) = 1 / cosh^2(x)
%     y = 1 ./ (cosh(x).^2);
% end
% 
% function ed = e_dot_placeholder(n)
%     % 若无法获得 \dot e，可置零；如需精确，请在外部计算并传入此函数中。
%     ed = zeros(n,1);
% end









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
DeltaU
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
t_eps = 1e-2;                 % 极小起跳时间
b_min = 1e-2;                 % 归一化时间下界，避免 log(0)
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









