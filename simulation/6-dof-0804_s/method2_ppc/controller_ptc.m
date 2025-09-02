

function  [dx,tau,alpha_sat,W_hat] = controller_ptc(t, x)
% t

n = 6;
q      = x(1:n);
dq     = x(n+1:2*n);
zeta1  = x(2*n+1:3*n);
zeta2  = x(3*n+1:4*n);
e1_int= x(4*n+1:5*n);
alpha_bar =x(5*n+1:6*n);   %36
c =x(6*n+1:6*n+2);%%38
d1=x(6*n+3:7*n+2); %
d2=x(7*n+3:8*n+2);

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
e1 = q - qd-zeta1;
e2 = dq - dqd-zeta2;


% 抗饱和项
T_zeta = 3;
gamma=0.4;
delta1 = zeros(n,1);
delta2 = zeros(n,1);
dzeta1 = zeros(n,1);
dzeta2 = zeros(n,1);
for i = 1:n
    delta1(i) = -12*zeta1(i);
    delta2(i) = -12*zeta2(i);
end

%误差转换
T_p =5;
p=0.6;
% kappa=1;
% dot_kappa=0;
[rho1, drho1] = performance_poly(t, T_p, pi/2, 0.02,p);
[rho2, drho2] = performance_poly(t, T_p, pi, 0.02,p);
% [kappa, dot_kappa] = shift(t, T_p,p);
kappa=1;dot_kappa=0;

% BLF函数自适应,如果设为1的话，无法在T_p=3s收敛
c_min=2;
c_max=6;
vartheta=0.1;    %自适应增益 γ_c

% Lyapunov 函数项
eps_den = 1e-4;   
%  disp(min(rho1.^2 - (kappa.*e1).^2));
den1 = max(rho1.^2 - (kappa.*e1).^2, eps_den);          % 限制下界
den2 = max(rho2.^2 - (kappa.*e2).^2, eps_den);          % 限制下界
V1 = (1/(2*c(1))) * sum( log( rho1.^2 ./ den1 ));
V2 = (1/(2*c(2))) * sum( log( rho2.^2 ./den2 ));
% % 自适应律
dot_c1 = -vartheta * V1;
dot_c2 = -vartheta * V2;

if c(1) <= c_min && dot_c1 < 0
    dot_c1 = 0;
elseif c(1) >= c_max && dot_c1 > 0
    dot_c1 = 0;
end

if c(2) <= c_min && dot_c2 < 0
    dot_c2 = 0;
elseif c(2) >= c_max && dot_c2 > 0
    dot_c2 = 0;
end
c(1) = min(max(c(1), c_min), c_max);
c(2) = min(max(c(2), c_min), c_max);

dot_c=[dot_c1;dot_c2];

% c(1)=2;
% c(2)=2;
% dot_c=[0;0];

% Control parameters
eta = 0.1;
iota_1 = 1;
iota_2 = 1;


k1 = diag([10, 10, 30, 8, 8, 10]*1);
k2 = diag([12, 10, 30, 10, 8, 8]*1);
kI = [12;12;30;12;12;12]*6;  
s1     = e2 + 8*e1;

% s1=8*e1+e2+theta;
sat_s1 = max(min(s1/0.1,1),-1);
hat_d1      = d1; 
ho=6;   w0=8;
dot_d1=-ho*sat_s1+ho/c(1)*(kappa*e1)./den1-w0*d1;
% dot_d1=-ho*sat_s1-w0*d1;
%—— 新的 alpha_raw—— - kI(i)*e1_int(i)
[h1,~] = rbfnn1([kappa*e1; kappa*e2], kappa*e1, rho1);
alpha_raw = zeros(n,1);
for i = 1:n
    alpha_raw(i) = dqd(i) - e2(i)+delta1(i)+dot_c(1)*e1(i)/(2*c(1)) ...
        -(iota_1 * dot_kappa^2 * e1(i)^3) /( den1(i) )...
        +e1(i)*(drho1/rho1)-k1(i) * e1(i)-hat_d1(i)- kI(i)*e1_int(i)-6*sat_s1(1)...
             - (c(1)*pi * e1(i)) / (2* eta * T_p )*(1+(n/2)^(eta/2)*((kappa*e1(i))^2 /  den1(i))^(eta/2));
%         - (c(1)*pi * e1(i)) / ( 2* eta * T_p )*(1+(n*kappa*e1(i)*den1(i))^(eta/2));

end
U_max1 = 100 * ones(6,1); 
% U_max1 = [60; 60; 40; 30; 20; 15]*0.85; 
alpha_sat=max( min(alpha_raw, U_max1), -U_max1 );
% alpha_sat=alpha_raw;
sat = @(x) x ./ max(1, abs(x));
mu_f = 0.4;  % 滤波时间常数
dalpha_bar = -( alpha_bar-alpha_sat) / mu_f;



s2=8*(q-alpha_sat)+(dq-dalpha_bar);
sat_s2 = max(min(s2/0.15,1),-1);

hat_d2      = d2;
dot_d2=-ho*sat_s1+ho/c(2)*(kappa*e2)./den2-w0*d2;
% dot_d2=-ho*sat_s2-w0*d2;
% hat_d2=-5*sat_s2(i);
h = zeros(n,1);
W_hat= zeros(n,1);
% RBFNN estimation and adaptive update
% [h,W_hat] = rbfnn([e1; e2], kappa*e2, rho2);
b = zeros(n,1);
for i = 1:n
    b(i) =dot_c(2)*e2(i)/(2*c(2)) - f(i)- h(i) + dalpha_bar(i)+delta2(i) ...
        -(iota_2 * dot_kappa^2 * e2(i)^3) /( den2(i))...
        + e2(i)*(drho2/rho2)-k2(i) * e2(i)-hat_d2(i)-6*sat_s1(1)...
        - (c(2)*pi * e2(i)) / (2* eta * T_p )*(1+(n/2)^(eta/2)*((kappa*e2(i))^2 /  den2(i))^(eta/2));
%       b(i) =dot_c(2)*e2(i)/(2*c(2)) - f(i)- h(i) - dalpha_bar(i)+delta2(i) ...
%         -(iota_2 * dot_kappa^2 * e2(i)^3) /( den2(i))...
%         + e2(i)*(drho2/rho2)-k2(i) * e2(i)-d2(i)...
%         - (c(2)*pi * e2(i)) / (2* eta * T_p )*(1+(n*kappa*e2(i)*den2(i))^(eta/2));
end

u=M*b;

% u_s =  tanh_saturation(u,U_max);
% U_max2 = 100 * ones(6,1); 
U_max2 = [60; 60; 40; 30; 20; 15]*1;  % 每个关节估算最大扭矩 [Nm]
u_s=max( min(u, U_max2), -U_max2 );
u_d=u-u_s;
tau = u_s;
ddq = f + M\tau + h;
%
% dzeta1 = delta1 + alpha_bar-alpha_raw + zeta2;
% dzeta2 =delta2 - M\u_d;

for i = 1:n
    dzeta1(i) = delta1(i)-12*zeta1(i) + alpha_bar(i) - alpha_raw(i) + zeta2(i);
    dzeta2(i) =delta2(i) -12*zeta2(i)- g(i,i)*u_d(i);
end
% dot_d1,dot_d2,d1,d2
if mod(round(t,2),0.5)==0      % 每 0.5 s 打一次
    fprintf('t=%.1f  d1=%+.2f  dot_d1=%+.3f  sat=%.2f\n',...
        t, d1(1), dot_d1(1), sat_s1(1));
end

dx = [dq; ddq; dzeta1;dzeta2;e1;dalpha_bar;dot_c;dot_d1;dot_d2];
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
function [kappa, dot_kappa] = shift(t, T_p, p)
    alpha = 6;
    beta_param = 4;

    % 初始化
    kappa = zeros(size(t));
    dot_kappa = zeros(size(t));

    % 完全Beta函数值
    B = beta(alpha, beta_param);

    % 归一化时间变量 beta(t)
    beta_t = (t ./ T_p).^p;

    % 拆分 t < T_p 和 t >= T_p
    idx1 = (t > 0) & (t < T_p);    % 处于上升阶段
    idx2 = ~idx1;        % 超出阶段

    % ---- 计算 \kappa(t) ----
    kappa(idx1) = betainc(beta_t(idx1), alpha, beta_param);  % 正则不完全Beta函数
    kappa(idx2) = 1;  % 饱和为1

    % ---- 计算 \dot{\kappa}(t) ----
    dbeta_dt = (p ./ T_p) .* (t ./ T_p).^(p - 1);  % \dot{\beta}(t)
    dot_kappa(idx1) = (1 ./ B) .* dbeta_dt(idx1) ...
        .* beta_t(idx1).^(alpha - 1) .* (1 - beta_t(idx1)).^(beta_param - 1);
    dot_kappa(idx2) = 0;
end



% function [kappa, dot_kappa] = shift(t, T_p)
%     p=0.98;
%     omega = [0 0.01 0.01 0.90 1.0];
%     n = 4;  % Bernstein基函数阶数
%     beta = min(t / T_p, 1);
% % 预分配存储
% B = zeros(n+1, length(t)); % Bernstein基函数 (i=0:4)
% B_dot = zeros(n+1, length(t)); % Bernstein基函数导数
% 
% % 计算Bernstein基函数及其导数
% for i = 0:n
%     B(i+1,:) = nchoosek(n,i) * beta.^i .* (1-beta).^(n-i);
%     B_dot(i+1,:) = nchoosek(n,i) * (i*beta.^(i-1).*(1-beta).^(n-i) - (n-i)*beta.^i.*(1-beta).^(n-i-1));
% end
% % 计算kappa和dot_kappa
% kappa = omega * B; % κ(t) = Σω_i B_{i,n}(β)
% dot_kappa = zeros(size(t));
% 
% for k = 1:length(t)
%     if t(k) >0 && t(k) < T_p
%         beta_dot = (p/T_p) * (t(k)/T_p)^(p-1);
%         dot_kappa(k) = omega * B_dot(:,k) * beta_dot;
%     else
%         kappa(k) = 1;
%         dot_kappa(k) = 0;
%     end
% end
% end




