

function  [dx,tau,alpha_sat,W_hat] = controller_ptc(t, x)
% t

n = 2;
q      = x(1:2);
dq     = x(3:4);
zeta1  = x(5:6);
zeta2  = x(7:8);
e1_int= x(9:10);
alpha_bar =x(11:12);
d1=x(13:14);
d2=x(15:16);

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
e1 = q - qd-zeta1;
e2 = dq - dqd-zeta2;


% 抗饱和项
T_zeta = 1.8;
gamma=0.4;
delta1 = zeros(n,1);
delta2 = zeros(n,1);
dzeta1 = zeros(n,1);
dzeta2 = zeros(n,1);
for i = 1:n
    delta1(i) = -12*zeta1(i)-pi/(4*T_zeta)*(abs(zeta1(i))^(1-gamma)+abs(zeta1(i))^(1+gamma));
    delta2(i) = -12*zeta2(i)-pi/(4*T_zeta)*(abs(zeta2(i))^(1-gamma)+abs(zeta2(i))^(1+gamma));
end

%误差转换
T_p = 3;
p=0.3;
% kappa=1;
% dot_kappa=0;
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
[rho1, drho1] = gppf(t, e1, DeltaU_now, d_est1, cfg1)
[rho2, drho2] = gppf(t, e2, DeltaU_now, d_est2, cfg2);
kappa=1;
dot_kappa=0;




% Lyapunov 函数项
eps_den = 1e-4;   
%  disp(min(rho1.^2 - (kappa.*e1).^2));
den1 = max(rho1.^2 - (kappa.*e1).^2, eps_den);          % 限制下界
den2 = max(rho2.^2 - (kappa.*e2).^2, eps_den);          % 限制下界


% Control parameters
eta = 0.1;
k1 = diag([12, 12]*3);
k2 = diag([6, 6]);
iota_1 = 1;
iota_2 = 1;
kI=[12, 12]*10;
c(1)=1;
c(2)=1;

s1  = e2 + 8*e1;

% s1=8*e1+e2+theta;
sat_s1 = max(min(s1/0.1,1),-1);
hat_d1      = d1; 
ho=6;   w0=8;
dot_d1=-ho*sat_s1+ho/(kappa*e1)./den1-w0*d1;

% [h1,~] = rbfnn1([kappa*e1; kappa*e2], kappa*e1, rho1);
alpha_raw = zeros(n,1);
for i = 1:n
    alpha_raw(i) = dqd(i) - e2(i)+delta1(i) ...
        -(iota_1 * dot_kappa^2 * e1(i)^3) /( den1(i) )...
        +e1(i)*(drho1(i)/rho1(i))-k1(i) * e1(i)-hat_d1(i)- kI(i)*e1_int(i)-6*sat_s1(1)...
             - (pi * e1(i)) / (2* eta * T_p )*(1+(n/2)^(eta/2)*((kappa*e1(i))^2 /  den1(i))^(eta/2));
%         - (pi * e1(i)) / ( 2* eta * T_p )*(1+(n*kappa*e1(i)*den1(i))^(eta/2));

end
U_max1 = [1000; 1000];  % 根据机械臂实际扭矩极限来选
alpha_sat=max( min(alpha_raw, U_max1), -U_max1 );
% alpha_sat=alpha_raw;
sat = @(x) x ./ max(1, abs(x));
mu_f = 0.2;  % 滤波时间常数
dalpha_bar = -( alpha_bar-alpha_sat) / mu_f;



s2=8*(q-alpha_sat)+(dq-dalpha_bar);
sat_s2 = max(min(s2/0.15,1),-1);

hat_d2      = d2;
dot_d2=-ho*sat_s1+ho/(kappa*e2)./den2-w0*d2;

h = zeros(2,1);
W_hat= zeros(2,1);
% RBFNN estimation and adaptive update
% [h,W_hat] = rbfnn([e1; e2], kappa*e2, rho2);
b = zeros(n,1);
for i = 1:2
    b(i) = - f(i)- h(i) + dalpha_bar(i)+delta2(i) ...
        -(iota_2 * dot_kappa^2 * e2(i)^3) /( den2(i))...
        + e2(i)*(drho2(i)/rho2(i))-k2(i) * e2(i)-hat_d2(i)-6*sat_s1(1)...
        - (pi * e2(i)) / (2* eta * T_p )*(1+(n/2)^(eta/2)*((kappa*e2(i))^2 /  den2(i))^(eta/2));
end

u=M*b;

% u_s =  tanh_saturation(u,U_max);
U_max2 = [10; 10]; 
u_s=max( min(u, U_max2), -U_max2 );
u_d=u-u_s;
tau = u_s;
ddq = f + M\tau + h;


for i = 1:n
    dzeta1(i) = delta1(i)-12*zeta1(i) + alpha_bar(i) - alpha_raw(i) + zeta2(i);
    dzeta2(i) =delta2(i) -12*zeta2(i)- g(i,i)*u_d(i);
end
% dot_d1,dot_d2,d1,d2
if mod(round(t,2),0.5)==0      % 每 0.5 s 打一次
    fprintf('t=%.1f  d1=%+.2f  dot_d1=%+.3f  sat=%.2f\n',...
        t, d1(1), dot_d1(1), sat_s1(1));
end

dx = [dq; ddq; dzeta1;dzeta2;e1;dalpha_bar;dot_d1;dot_d2];
end

function [rho, drho] = performance_poly(t, T_p, rho0, a,p)
if t>0
    beta0 = (t ./ T_p).^p;
    beta0 = min(max(beta0, 0), 1);

    % 在 [0, T_p] 内计算四次多项式值
    h    = 1 - 4*beta0.^3 + 3*beta0.^4;
    rho  = a + (rho0 - a) .* h;

    % 对应导数
    if nargout > 1
        % d/dβ0 [1 - 4β0^3 + 3β0^4] = -12β0^2 + 12β0^3
        dhdB0 = -12*beta0.^2 + 12*beta0.^3;
        % dβ0/dt = p/T_p * (t/T_p)^(p-1) = (p/T_p)*beta0.^(1 - 1/p)
        dB0dt = (p ./ T_p) .* ( (t ./ T_p).^(p - 1) );
        % 所以 drho
        drho = (rho0 - a) .* dhdB0 .* dB0dt;
        % t > T_p 时，强制 drho=0
        drho(t >= T_p) = 0;
    end
else
    rho=rho0;
    drho=0;
end
end








