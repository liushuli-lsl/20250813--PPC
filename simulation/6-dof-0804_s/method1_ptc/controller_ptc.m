

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
T_zeta = 1.5;
gamma=0.4;
delta1 = zeros(n,1);
delta2 = zeros(n,1);
dzeta1 = zeros(n,1);
dzeta2 = zeros(n,1);
for i = 1:n
    delta1(i) = -8*zeta1(i)-pi/(4*T_zeta)*(abs(zeta1(i))^(1-gamma)+abs(zeta1(i))^(1+gamma));
    delta2(i) = -8*zeta2(i)-pi/(4*T_zeta)*(abs(zeta2(i))^(1-gamma)+abs(zeta2(i))^(1+gamma));
end

%误差转换
T_p = 3.5;

c(1)=1;
c(2)=1;
dot_c=[0;0];

% Control parameters
eta = 0.2;

k1 = diag([10, 10, 10, 8, 8, 10]*30);
k2 = diag([12, 10, 12, 10, 8, 8]*2);
kI = [12;12;12;12;12;12]*2;  
s1     = e2 + 8*e1;

% s1=8*e1+e2+theta;
sat_s1 = max(min(s1/0.1,1),-1);
hat_d1      = d1; 
ho=6;   w0=8;
dot_d1=-ho*sat_s1-w0*d1;
scale = max((T_p - t) / T_p, 0.1);  % 最小缩放限制，防止分母趋零
% 虚拟控制律 alpha(t)
eps = 1e-3;
alpha_raw = -k1 * e1 ./ ((scale).^eta)+delta1- kI(i)*e1_int(i) - hat_d1;
U_max1 = 500 * ones(n,1);
alpha_sat = max(min(alpha_raw, U_max1), -U_max1);

% 滤波器：alpha_bar(t)
mu_f = 0.5;
dalpha_bar = -(alpha_bar - alpha_sat) / mu_f;

% 扰动估计器 d2
hat_d2 = d2;
dot_d2 = -ho * sat_s1 - w0 * d2;
k1b = diag([1, 2, 1, 1.5, 1.2, 1]*40);   % 用于 b 控制律中的辅助位置反馈
% 控制律 tau(t)
b = -f - k2 * e2 ./ ((scale).^eta) - hat_d2 + dalpha_bar + delta2-k1b * e1;
u = M * b;
h = zeros(n,1);
W_hat = zeros(n,1);
% u_s =  tanh_saturation(u,U_max);
% U_max2 = 100 * ones(6,1); 
U_max2 = [60; 60; 40; 30; 20; 15]*1;  % 每个关节估算最大扭矩 [Nm]
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

dx = [dq; ddq; dzeta1;dzeta2;e1;dalpha_bar;dot_c;dot_d1;dot_d2];
end

