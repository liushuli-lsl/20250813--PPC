

function  [dx,tau,alpha_sat] = controller_ptc(t, x)
% t

n = 2;
q      = x(1:2);
dq     = x(3:4);
zeta  = x(5:6);
z1_int  = x(7:8);
alpha_bar= x(9:10);
d1 =x(11:12);
d2=x(13:14);


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
T_p = 3;
p=0.3;
% [rho1, drho1] = performance_poly(t, T_p, 2, 0.02,p);
% [rho2, drho2] = performance_poly(t, T_p, 4, 0.01,p);


DeltaU_now = zeros(n,1);    % 当前步尚未知饱和残差，先占位
d_est1     = abs(d1);       % 你已有的 d1 自适应估计
d_est2     = abs(d2);

% 配置不同通道的 G-PPF 参数（可按需微调）
cfg1 = struct('id',1,'Tp',T_p,'p',p,'a',0.02, ...   % 位置误差通道下界 a1
    'sigma0',1.0,'sigma_min',0.6,'sigma_max',1.5, ...
    'iota',3.0,'Sigma_max',0.5, ...
    'k_u',0.7,'k_d',0.4,'k_e',0.2, ...
    'use_lpf',true,'tau_u',0.6,'tau_d',0.8,'tau_e',0.8);

cfg2 = struct('id',2,'Tp',T_p,'p',p,'a',0.01, ...   % 速度误差通道下界 a2
    'sigma0',1.0,'sigma_min',0.6,'sigma_max',1.5, ...
    'iota',3.0,'Sigma_max',0.5, ...
    'k_u',0.7,'k_d',0.4,'k_e',0.2, ...
    'use_lpf',true,'tau_u',0.6,'tau_d',0.8,'tau_e',0.8);

% 计算 G-PPF 及其导数（用于 BLF 精确补偿）
[rho1, drho1] = gppf(t, z1, DeltaU_now, d_est1, cfg1);
[rho2, drho2] = gppf(t, z2, DeltaU_now, d_est2, cfg2);

% Lyapunov 函数项
eps_den = 1e-4;
%  disp(min(rho1.^2 - (z1).^2));
den1 = zeros(n,1);
den2 = zeros(n,1);
for i=1:n
    den1(i) = max(rho1(i).^2 - (z1(i)).^2, eps_den);          % 限制下界den
    den2(i) = max(rho2(i).^2 - (z2(i)).^2, eps_den);          % 限制下界
end

% Control parameters
eta = 0.1;
k1 = diag([12, 12]*3);
k2 = diag([6, 6]);
kI=[12, 12]*10;

hat_d1      = d1;
ho=6;   w0=8;

dot_d1= zeros(n,1);
dot_d2= zeros(n,1);
s1= zeros(n,1);
sat_s1= zeros(n,1);
alpha_raw = zeros(n,1);
for i = 1:n
    s1(i)  = z2(i) + 8*z1(i);
    sat_s1(i) = max(min(s1(i)/0.1,1),-1);
    dot_d1(i)=-ho*sat_s1(i)+ho/(z1(i))./den1(i)-w0*d1(i);
    alpha_raw(i) = dqd(i) - z2(i) ...
        -( z1(i)^3) /( den1(i) )...
        +z1(i)*(drho1(i)/rho1(i))-k1(i) * z1(i)-hat_d1(i)- kI(i)*z1_int(i)-6*sat_s1(1)...
        - (pi * z1(i)) / (2* eta * T_p )*(1+(n/2)^(eta/2)*((z1(i))^2 /  den1(i))^(eta/2));
    %         - (pi * z1(i)) / ( 2* eta * T_p )*(1+(n*z1(i)*den1(i))^(eta/2));

end
% U_max1 = [1000; 1000];  % 根据机械臂实际扭矩极限来选
alpha_sat=0;
mu_f = 0.2;  % 滤波时间常数
dalpha_bar = -( alpha_bar-alpha_raw) / mu_f;


% 抗饱和项
T_zeta = 1.8;
gamma=0.4;
sigpow = @(x,p) sign(x).*abs(x).^p;
delta = zeros(n,1);
dzeta = zeros(n,1);
for i = 1:n
%     delta(i) = -12*zeta(i)-pi/(4*T_zeta)*(abs(zeta(i))^(1-gamma)+abs(zeta(i))^(1+gamma));
delta  =-12*zeta(i) -(pi/(2*gamma*T_zeta)) * ( sigpow(zeta,1-gamma) + sigpow(zeta,1+gamma) );
end


hat_d2      = d2;
h= zeros(n,1);
b = zeros(n,1);
for i = 1:2
    dot_d2(i)=-ho*sat_s1(i)+ho/(z2(i))./den2(i)-w0*d2(i);
    b(i) = - f(i)- h(i) + dalpha_bar(i)+delta(i) ...
        -( z2(i)^3) /( den2(i))...
        + z2(i)*(drho2(i)/rho2(i))-k2(i) * z2(i)-hat_d2(i)-6*sat_s1(1)...
        - (pi * z2(i)) / (2* eta * T_p )*(1+(n/2)^(eta/2)*((z2(i))^2 /  den2(i))^(eta/2));
end

u=M*b;

% u_s =  tanh_saturation(u,U_max);
U_max2 = [10; 10];
tau=max( min(u, U_max2), -U_max2 );
u_d=u-tau;
ddq = f + M\tau + h;

% 回灌 Δu/d 到 G-PPF 的 LPF（只更新内部状态，不重算本步 rho）
cfg1.feed_only = true;  gppf(t, z1, u_d, d_est1, cfg1);
cfg2.feed_only = true;  gppf(t, z2, u_d, d_est2, cfg2);


for i = 1:n
    dzeta(i) =delta(i) -12*zeta(i)- g(i,i)*u_d(i);
end
if mod(round(t,2),0.5)==0      % 每 0.5 s 打一次
    fprintf('t=%.1f  d1=%+.2f  dot_d1=%+.3f  sat=%.2f\n',...
        t, d1(1), dot_d1(1), sat_s1(1));
end


dx = [dq; ddq; dzeta;z1;dalpha_bar;dot_d1;dot_d2];
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








