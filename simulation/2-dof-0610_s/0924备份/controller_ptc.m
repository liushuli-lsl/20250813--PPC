

function  [dx,tau,rho1,w] = controller_ptc(t, x)
% t

n = 2;
q      = x(1:2);
dq     = x(3:4);
zeta  = x(5:6);
z1_int  = x(7:8);
alpha_bar= x(9:10);
% d2=x(13:14);
theta=x(11:12);
omegahat=x(13:14);

qd = [0.1*sin(0.5*t) + cos(0.5*t);0.1*sin(t) + cos(t)];
dqd = [0.05*cos(0.5*t)-0.5*sin(0.5*t); 0.1*cos(t)-sin(t)];
ddqd = [-0.025*sin(0.5*t)-0.25*cos(0.5*t);
    -0.1*sin(   t)-   cos(   t)];

% 系統方程
% Dynamics
[M, C, G] = two_link_dynamics(q, dq);
f = -M \ (C*dq + G);
g = M \ eye(2);

% Errors
z1 = q - qd;
z2 = dq - dqd-zeta;


%误差转换
T_p = 2;
p=0.5;
% [rho1, drho1] = performance_poly(t, T_p, 2, 0.01,p);
% [rho2, drho2] = performance_poly(t, T_p, 4, 0.01,p);
a=0.01;

DeltaU_now = zeros(n,1);    % 当前步尚未知饱和残差，先占位
d_est1     = abs(omegahat);       % 你已有的 theta 自适应估计
% d_est2     = abs(d2);

% 配置不同通道的 G-PPF 参数（可按需微调）
cfg1 = struct('id',1,'Tp',T_p,'p',p,'a',a, ...   % 位置误差通道下界 a1
    'sigma0',1.0,'sigma_min',1.2,'sigma_max',2.2, ...
    'iota',3.0,'Sigma_max',0.5, ...
    'k_u',1.5,'k_d',2.5,'k_e',0.4, ...
    'use_lpf',true,'tau_u',0.2,'tau_d',0.2,'tau_e',0.4);


% Sigma_max  = 1.00;     % 最大放宽倍数：rho = a*(1 + g*Sigma)（总幅 ~ a*(1+Sigma_max)）
% k_d        = 2.50;     % 扰动权重
% tau_d_post = 0.15;     % 扰动低通
% k_e        = 1.20;     % 误差权重
% tau_e      = 0.07;     % 误差低通

% 计算 G-PPF 及其导数（用于 BLF 精确补偿）
[rho1, drho1] = gppf(t,n, z1, [], d_est1, cfg1);

eps0 = 1e-6;
eta  = 0.1;
% Control parameters
k1 = diag([12, 12])*1;
k2 = diag([6, 6])*0.1;
Lambda=diag([5, 5]);
beta = 0.2; 
omega_d = 8;                 % 标量或逐通道统一标量
k_s  = 0; 

rho1_safe = max(rho1,1e-6);
r1 = min(0.98, max(-0.98, z1./rho1_safe));
% 相对下界：den >= 0.02*rho^2
den1 = max(rho1.^2 - z1.^2, a*rho1.^2 + eps0);
[P1,~,V1] = blf_terms(z1, rho1);
K1   = (2*pi/(eta*T_p)) * ( max(V1,eps0).^( - eta/2) + (n^(-eta/2)) * max(V1,eps0).^(eta/2) );         % \mathcal{K}_{1}
% 虚拟控制器
% ===== 调用扰动模块 =====
mode  = 'both';                 % 'none' | 'matched' | 'unmatched' | 'both'
level = [20;20];                  % 扰动幅值
band  = [8; 8];                  % 基波频率
seed  = 0;                       % 随机种子
[w, nu, d_true,tau_d] = disturbance(t, q, dq, mode, level, band, seed);
Phi =(den1).*z1./rho1.^4;
alpha = dqd - zeta ...
      + (r1.^3) .* drho1  - k1* Phi ;    % (z1^3/rho1^3)*drho1        - K1.* Phi1 - k1* Phi1 ..

      U_alpha = [100;100];   % 机械臂虚拟输入允许范围
% alpha_sat = max(min(alpha, U_alpha), -U_alpha);                  % 0<beta<1
dalpha_bar = -(alpha_bar - alpha) / beta;


%% --------- 7) 动态抗饱和补偿 \zeta ----------
T_zeta = 1.8;  gamma = 0.4;  mu = 12;        % 0<gamma<1, mu2>0
sigpow = @(x,p) sign(x).*abs(x).^p;
delta  = -(pi/(2*gamma*T_zeta)) * ( sigpow(zeta,1-gamma) + sigpow(zeta,1+gamma) );

%% --------- 8) Step-2：实际控制扭矩 u（按文稿式） ----------
K2   = (pi/(2*eta*T_p)) *  ( 2^(eta/2) * sign(z2) .* abs(z2).^(1-eta) ...
              + (n/2)^(eta/2) * sign(z2) .* abs(z2).^(1+eta) );            % \mathcal{K}_{2}


% 未建模项估计（若无 RBFNN，则置零）
% omegahat = zeros(n,1);                                         % sgn 增益
u = C*dq + G + M*( dalpha_bar + delta - mu*zeta ...
        - P1.* z1  - K2  - omegahat  ...
        - k_s * sign(z2)  ...
        - k2 * z2);

% if mod(round(t,2),1)==0      % 每 0.5 s 打一次
%     fprintf('t=%.1f  alpha=%+.2f   u=%+.2f\n',...
%         t, alpha, u);
% end
% u_s =  tanh_saturation(u,U_max);
U_max2 = [10; 10];
tau=max( min(u, U_max2), -U_max2 );
u_d=u-tau;



  
ddq = f + g*(tau)+w ;

%     一阶跟踪微分器（避免直接微分dq）
    theta_dot   = -omega_d*theta + omega_d*dq;
    x2hat_dot   =  omega_d*(dq - theta);
    omegahat_dot = -Lambda*omegahat + Lambda*( x2hat_dot - f - g*tau );
 
% 回灌 Δu/d 到 G-PPF 的 LPF（只更新内部状态，不重算本步 rho）
u_d_abs = min(max(abs(u_d), 0), 15)        % |Δu| 逐通道，上限按实际设
d_abs   = min(max(abs(omegahat), 0), 200)   % |扰动估计|；验证时也可用 abs(d_true)
cfg1.feed_only = true;  gppf(t, n,z1,u_d_abs, d_abs, cfg1);

% ddq = f + M\tau +omegahat+nu+w;
% 
% %     一阶跟踪微分器（避免直接微分dq）
%     theta_dot   = -omega_d*theta + omega_d*dq;
%     x2hat_dot   =  omega_d*(dq - theta);
%     omegahat_dot = -Lambda*omegahat + Lambda*( x2hat_dot - f - g*tau );
%  
% % 回灌 Δu/d 到 G-PPF 的 LPF（只更新内部状态，不重算本步 rho）
% u_d_abs = min(max(abs(u_d), 0), 15)        % |Δu| 逐通道，上限按实际设
% d_abs   = min(max(abs(omegahat), 0), 20)   % |扰动估计|；验证时也可用 abs(d_true)
% cfg1.feed_only = true;  gppf(t, n,z1,u_d_abs, d_abs, cfg1);

dzeta = delta - mu*zeta - g*u_d;
dx = [dq; ddq; dzeta;z1;dalpha_bar;theta_dot;omegahat_dot];
end


%% ======= 局部小工具：BLF 基元与 Psi 函数 =======
function [P,Q,V] = blf_terms(z, rho)
epsd = 1e-6;
a=0.002;
% den  = max(rho.^2 - z.^2, epsd);
den = max(rho.^2 - z.^2, a*rho.^2 + epsd);
P = (rho.^4) ./ max(den.^2, epsd);
Q = (z.^4)   ./ max(den.^2, epsd);
V = 0.5 * (rho.^2 .* z.^2) ./ den;
end

