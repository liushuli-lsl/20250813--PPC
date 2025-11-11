

function  [dx,tau,rho1,w,u] = controller_ptc(t, x)
% t

n = 6;
q      = x(1:n);
dq     = x(n+1:2*n);
zeta  = x(2*n+1:3*n);
zeta2  = x(3*n+1:4*n);
alpha_bar= x(4*n+1:5*n);
theta =x(5*n+1:6*n);   %36
omegahat =x(6*n+1:7*n);%%38
% =x(6*n+3:7*n+2); %


% 位置轨迹 (组合正弦和余弦)
qd = [0.1*sin(0.5*t) + cos(0.5*t);    % 关节1（低频主导）
      0.1*sin(t) + cos(t);            % 关节2（中频）
      0.2*sin(1.5*t) + 0.8*cos(t);    % 关节3（混合频率）
      0.3*sin(2*t) + 0.7*cos(0.5*t);  % 关节4（高频+低频）
      0.1*sin(0.3*t) + 0.9*cos(0.2*t);% 关节5（超低频）
      0.4*sin(t) + 0.6*cos(2*t)];     % 关节6（交叉频率）

% 速度轨迹（解析求导）
dqd = [0.05*cos(0.5*t) - 0.5*sin(0.5*t);   % 关节1
       0.1*cos(t) - sin(t);                % 关节2
       0.3*cos(1.5*t) - 0.8*sin(t);        % 关节3
       0.6*cos(2*t) - 0.35*sin(0.5*t);     % 关节4
       0.03*cos(0.3*t) - 0.18*sin(0.2*t);  % 关节5
       0.4*cos(t) - 1.2*sin(2*t)];         % 关节6

% 系統方程
% Dynamics
[M, C, G] = six_link_dynamics(q, dq);
f = -M \ (C + G);
g = M \ eye(n);

% Errors
z1 = q - qd;
z2 = dq - dqd-zeta;


%误差转换
T_p = 3; a=0.02;
p=0.6;
% DeltaU_now = zeros(n,1);    % 当前步尚未知饱和残差，先占位
d_est1     = abs(omegahat);       % 你已有的 theta 自适应估计
% d_est2     = abs(d2);

% 配置不同通道的 G-PPF 参数（可按需微调）
cfg = struct('id',1,'Tp',T_p,'p',p,'a',a, ...   % 位置误差通道下界 a1
    'sigma0',1.0,'sigma_min',1.2,'sigma_max',2, ...
    'iota',2,'Sigma_max',20, ...
    'k_u',0.2,'k_d',0.15,'k_e',0.4, ...
    'use_lpf',true,'tau_u',0.05,'tau_d',0.1,'tau_e',0.08);

persistent Ud_prev d_prev init
if isempty(init), Ud_prev = zeros(n,1); d_prev = zeros(n,1); init = true; end

cfg.feed_only = true;
gppf(t, n, z1, Ud_prev, d_prev, cfg);  % 只更新LPF，不取rho
cfg = rmfield(cfg, 'feed_only');       % 避免污染后续
[rho1, drho1] = gppf(t,n, z1, Ud_prev, d_prev,cfg);
% [rho1, drho1] = performance_poly(t, T_p, pi/2, 0.02,p);


% Control parameters
eps0 = 1e-6;
eta = 0.1;

k1 = diag([10, 10, 10, 8, 8, 10]*6);
k2 = diag([12, 10, 12, 10, 8, 8]*1);
kI = [12;12;12;12;12;12]*20;  


Lambda=diag([5,5,5,5,5, 5])*12;
beta = 0.4; 
omega_d = 30;                 % 标量或逐通道统一标量
k_s  = 0;
rho1_safe = max(rho1,1e-6);
r1 = min(0.98, max(-0.98, z1./rho1_safe));
% 相对下界：den >= 0.02*rho^2
den1 = max(rho1.^2 - z1.^2, a*rho1.^2 + eps0);
[P1,~,V1] = blf_terms(z1, rho1);
K1   = (pi/(2*eta*T_p)) * ( max(V1,eps0).^( - eta/2) + (n^(-eta/2)) * max(V1,eps0).^(eta/2) ).*rho1.^2;         % \mathcal{K}_{1}
% 虚拟控制器*

Phi =(den1).*z1./rho1.^4;
alpha = dqd - zeta ...
      + (r1.^3) .* drho1  - K1.* Phi - k1* Phi;     % (z1^3/rho1^3)*drho1        - K1.* Phi1 - k1* Phi1 ..
% alpha_sat = max(min(alpha, U_alpha), -U_alpha);                  % 0<beta<1
dalpha_bar = -(alpha_bar - alpha) / beta;

% 抗饱和项
T_zeta = 1.8;
gamma=0.4; mu = 14;
sigpow = @(x,p) sign(x).*abs(x).^p;
delta  = -(pi/(2*gamma*T_zeta)) * ( sigpow(zeta,1-gamma) + sigpow(zeta,1+gamma) );

%% --------- 8) Step-2：实际控制扭矩 u（按文稿式） ----------
K2   = (pi/(2*eta*T_p)) *  ( 2^(eta/2) * sign(z2) .* abs(z2).^(1-eta) ...
              + (n/2)^(eta/2) * sign(z2) .* abs(z2).^(1+eta) );            % \mathcal{K}_{2}
                                       % sgn 增益
u = C.*dq + G + M*( dalpha_bar + delta - mu*zeta ...
        - P1.* z1  - K2 ...  
        - k2 * z2);
% U_max2 = [60; 60; 40; 30; 30; 10]*3;  % 每个关节估算最大扭矩 [Nm]
U_max2 =[100; 100; 40; 30; 20; 15];
tau=max( min(u, U_max2), -U_max2 );
u_d=u-tau;


% % ===== 调用扰动模块 =====
mode  = 'both';                 % 'none' | 'matched' | 'unmatched' | 'both'
level = [60;60;60;60;60;60]*0.125; % 扰动幅值
band  = [4; 4];                  % 基波频率
seed  = 0;                       % 随机种子
[w, nu, d_true,tau_d] = disturbance(t, q, dq, mode, level, band, seed);
%   w=zeros(n,1);
ddq = f + g*(tau) +w;

%     一阶跟踪微分器（避免直接微分dq）
    theta_dot   = -omega_d*theta + omega_d*dq;
    x2hat_dot   =  omega_d*(dq - theta);
    omegahat_dot = -Lambda*omegahat + Lambda*( x2hat_dot - f - g*tau );

 Ud_prev = abs(u_d);           % 本步的 |ΔU|
d_prev  = abs(omegahat);          % 或 abs(d_true)

dzeta = delta - mu*zeta - g*u_d;
dx = [dq; ddq; dzeta;z1;dalpha_bar;theta_dot;omegahat_dot];
% dx = [dq; ddq; dzeta1;dzeta2;e1;dalpha_bar;dot_c;dot_d1;dot_d2];
end


%% ======= 局部小工具：BLF 基元与 Psi 函数 =======
function [P,Q,V] = blf_terms(z, rho)
epsd = 1e-6;
a=0.02;
% den  = max(rho.^2 - z.^2, epsd);
den = max(rho.^2 - z.^2, a*rho.^2 + epsd);
P = (rho.^4) ./ max(den.^2, epsd);
Q = (z.^4)   ./ max(den.^2, epsd);
V = 0.5 * (rho.^2 .* z.^2) ./ den;
end

% function [rho, drho] = performance_poly(t, T_p, rho0, a,p)
% if t>0
%     beta0 = (t ./ T_p).^p;
%     beta0 = min(max(beta0, 0), 1);
% 
%     % 在 [0, T_p] 内计算四次多项式值
%     h    = 1 - 4*beta0.^3 + 3*beta0.^4;
%     rho  = a + (rho0 - a) .* h;
% 
%     % 对应导数
%     if nargout > 1
%         % d/dβ0 [1 - 4β0^3 + 3β0^4] = -12β0^2 + 12β0^3
%         dhdB0 = -12*beta0.^2 + 12*beta0.^3;
%         % dβ0/dt = p/T_p * (t/T_p)^(p-1) = (p/T_p)*beta0.^(1 - 1/p)
%         dB0dt = (p ./ T_p) .* ( (t ./ T_p).^(p - 1) );
%         % 所以 drho
%         drho = (rho0 - a) .* dhdB0 .* dB0dt;
%         % t > T_p 时，强制 drho=0
%         drho(t >= T_p) = 0;
%     end
% else
%     rho=rho0;
%     drho=0;
% end
% end
% function [kappa, dot_kappa] = shift(t, T_p, p)
%     alpha = 6;
%     beta_param = 4;
% 
%     % 初始化
%     kappa = zeros(size(t));
%     dot_kappa = zeros(size(t));
% 
%     % 完全Beta函数值
%     B = beta(alpha, beta_param);
% 
%     % 归一化时间变量 beta(t)
%     beta_t = (t ./ T_p).^p;
% 
%     % 拆分 t < T_p 和 t >= T_p
%     idx1 = (t > 0) & (t < T_p);    % 处于上升阶段
%     idx2 = ~idx1;        % 超出阶段
% 
%     % ---- 计算 \kappa(t) ----
%     kappa(idx1) = betainc(beta_t(idx1), alpha, beta_param);  % 正则不完全Beta函数
%     kappa(idx2) = 1;  % 饱和为1
% 
%     % ---- 计算 \dot{\kappa}(t) ----
%     dbeta_dt = (p ./ T_p) .* (t ./ T_p).^(p - 1);  % \dot{\beta}(t)
%     dot_kappa(idx1) = (1 ./ B) .* dbeta_dt(idx1) ...
%         .* beta_t(idx1).^(alpha - 1) .* (1 - beta_t(idx1)).^(beta_param - 1);
%     dot_kappa(idx2) = 0;
% end





