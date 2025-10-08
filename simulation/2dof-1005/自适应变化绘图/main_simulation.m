clearvars;          % 清除工作区变量
clear functions;    % 清除所有函数的 persistent 缓存
close all;          % 关闭所有图窗
clc;                % 清空命令窗口
clear all;


%% 1) 仿真设置
n=2;  %自由度
T_total = 20; 
N       = 20000;            % 点数
tspan   = linspace(0, T_total, N);
dt      = tspan(2) - tspan(1);

%% 2) 初始状态
q00 = [
          pi,     pi;         % 偏离参考较大
    pi/2,     pi/2;        % 中等伸展                     T=3.5,kI=[12, 6]*10;
0,    pi/2;        % 十字形对称
 pi/4,     pi/4;        % 十字形对称
      -pi/2,    -pi/2;        % 十字形对称
       -3*pi/4,     -3*pi/4;      % 偏离参考较大
          0,        0;          % 折叠（近奇异）
];
q0 = q00(6, :)';   % q0 = [0; 0]  2,3,4,5
dq0 = [0;0];                  %初始誤差不能爲0，否則會產生歧義
zeta   = zeros(n,1);
e1_int  = zeros(n,1);
alpha_bar = zeros(n,1);
theta = zeros(n,1);
omegahat = zeros(n,1);
x0  = [q0; dq0;zeta;e1_int;alpha_bar;theta;omegahat];


%% 3) 预分配存储

x       = zeros(7*n, N);
q_use   = zeros(N,n);
dq_use  = zeros(N,n);
qd_mat  = zeros(N,n);
dqd_mat = zeros(N,n);
tau_mat = zeros(N,n);
rho_mat = inf(N, n);      % 全部初始化为 +Inf
alpha_mat = zeros(N,n);
zeta_mat= zeros(N,n);
omegahat_mat = zeros(N,n);
d_mat = zeros(N,n);
u_mat = zeros(N,n);
rho_p_mat=zeros(N,n);
rho_m_mat=zeros(N,n);

x(:,1) = x0;
q_use(1,:)  = q0.';
dq_use(1,:) = dq0.';

%% 4) 固定步长 RK4 主循环
for k = 1:N-1
    t = tspan(k)
    % 解算一次 RK4
  [k1, tau1,rho1,d_true,u1] = controller_ptc(t,           x(:,k));
[k2, tau2,~,~,~] = controller_ptc(t+dt/2,      x(:,k)+dt/2*k1);
[k3, tau3,~,~,~] = controller_ptc(t+dt/2,      x(:,k)+dt/2*k2);
[k4, tau4,~,~,~] = controller_ptc(t+dt,        x(:,k)+dt*k3);

    x(:,k+1) = x(:,k) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;

    % 存储关节、参考、控制
    q   = x(1:2,k+1);
    dq  = x(3:4,k+1);
    zeta= x(5:6,k+1);
     alpha= x(9:10,k+1);
omegahat=x(13:14,k+1);
% qd = [0.1*sin(0.5*t) + cos(0.5*t);0.1*sin(t) + cos(t)];
% dqd = [0.05*cos(0.5*t)-0.5*sin(0.5*t); 0.1*cos(t)-sin(t)];
% ===== 提高参考轨迹频率：用统一缩放系数 s (>1) =====
s = 3.0;                 % 频率放大倍数：2=加倍，3=三倍... 可自行调整
w1 = 1.5*s;              % 关节1的新角频率
w2 = 1.5*s;              % 关节2的新角频率

qd   = [ 0.1*sin(w1*t) + cos(w1*t);
         0.1*sin(w2*t) + cos(w2*t) ];

% 一阶导（速度）：d/dt[sin(wt)]=w*cos(wt), d/dt[cos(wt)]=-w*sin(wt)
dqd  = [ w1*( 0.1*cos(w1*t) - sin(w1*t) );
         w2*( 0.1*cos(w2*t) - sin(w2*t) ) ];

    q_use(k+1,:)   = q.';
    dq_use(k+1,:)  = dq.';
    qd_mat(k,:)  = qd.';
    dqd_mat(k,:) = dqd.';

    tau_mat(k+1,:) = tau1.';
     alpha_mat(k+1,:) = alpha.';
     zeta_mat(k+1,:) = zeta.';
 rho_mat(k+1,:) = rho1.';
   omegahat_mat(k+1,:) =  omegahat.';
     d_mat(k+1,:) =  d_true.';
       u_mat(k+1,:) = u1.';
% rho_p_mat(k+1,:) = rho_p.';
% rho_m_mat(k+1,:) = rho_m.';
       e=zeros(n,1);  
uu = zeros(n,1);    % 当前步尚未知饱和残差，先占位
d     = zeros(n,1);       % 你已有的 theta 自适应估计
% 配置不同通道的 G-PPF 参数（可按需微调）
cfg2 = struct('id',1,'Tp',3,'p',0.3,'a',0.01, ...   % 位置误差通道下界 a1
    'sigma0',1.0,'sigma_min',1.2,'sigma_max',2, ...
    'iota',3,'Sigma_max',20, ...
    'k_u',1.8,'k_d',0.3,'k_e',0.8, ...
    'use_lpf',true,'tau_u',0.05,'tau_d',0.1,'tau_e',0.08);
% % 计算 G-PPF 及其导数（用于 BLF 精确补偿）
% [rho_p2, rho_m2,rho_dot_p2,rho_dot_m2] =   gppf2(t,n, e, uu, d,cfg2);

end

e_q  = q_use   - qd_mat;
e_dq = dq_use  - dqd_mat;



% [rho1] = arrayfun(@(tt) performance_poly1(tt),tspan);
% [rho2] = arrayfun(@(tt) performance_poly2(tt),tspan);

save('method1.mat', 'tspan', 'e_q', 'e_dq',  'qd_mat', 'q_use','dqd_mat', 'dq_use','rho_mat','tau_mat','u_mat',"alpha_mat",'omegahat_mat');
figure;
for i = 1:n
    subplot(2,1,i);
    plot(tspan, e_q(:,i), 'b', 'LineWidth', 1.5); hold on;
    plot(tspan, rho_mat(:,i), 'k-', 'LineWidth',1.5);hold on;
      plot(tspan,  -rho_mat(:,i), 'k-', 'LineWidth',1.5);
    yline(0, 'k--');
    xline(3, 'r--', 'LineWidth', 1.2);
    title(['Tracking Error e_' num2str(i)]);
    xlabel('Time (s)'); ylabel('e_i (rad)');
    legend('Error', 'Zero Line', 'T_p');
    ylim([-0.5,0.5]);
end


% Plot error curves with T_p marker
figure;
for i = 1:n
    subplot(2,1,i);
    plot(tspan, e_dq(:,i), 'b', 'LineWidth', 1.5); hold on;
%       plot(tspan, rho2, 'k--', 'LineWidth',1.5);hold on;
%       plot(tspan, -rho2, 'k--', 'LineWidth',1.5);
    yline(0, 'k--');
    xline(3, 'r--', 'LineWidth', 1.2);
    title(['Tracking Error e_' num2str(i)]);
    xlabel('Time (s)'); ylabel('e_i (rad)');
    legend('Error', 'Zero Line', 'T_p');
    ylim([-0.5,0.5]);
end





% 
% 
% % Plot joint position vs reference
figure;
for i = 1:n
    subplot(2,1,i);
    plot(tspan, q_use(:,i), 'b', 'LineWidth', 1.5); hold on;
    plot(tspan, qd_mat(:,i), 'r--', 'LineWidth', 1.2);
    title(['Joint q_' num2str(i)]);
    legend('Actual', 'Reference');
    xlabel('Time (s)'); ylabel('Position (rad)');
end




% Plot joint velocity vs reference
% figure;
% for i = 1:n
%     subplot(2,1,i);
%     plot(tspan, dq_use(:,i), 'b', 'LineWidth', 1.5); hold on;
%     plot(tspan, dqd_mat(:,i), 'r--', 'LineWidth', 1.2);
%     title(['Joint dq_' num2str(i)]);
%     legend('Actual', 'Reference');
%     xlabel('Time (s)'); ylabel('Velocity (rad/s)');
% end
% % 
% % Plot control torques
figure;
for i = 1:n
    subplot(2,1,i);
      plot(tspan, u_mat(:,i), 'r-', 'LineWidth', 1.5); hold on;
  
    plot(tspan, tau_mat(:,i), 'b--', 'LineWidth', 1.5);
    title(['Joint \tau_' num2str(i)]);
    xlabel('Time (s)'); ylabel('Torque (Nm)');
    ylim([-10,10]);
end
% figure;
% for i = 1:n
%     subplot(2,1,i);0
%     plot(tspan, alpha_mat(:,i), 'b', 'LineWidth', 1.5);
%     title(['Joint \alpha_' num2str(i)]);
%     xlabel('Time (s)'); ylabel('Torque (Nm)');
% end
% figure;
% for i = 1:n
%     subplot(2,1,i);
%     plot(tspan, omegahat_mat(:,i), 'k', 'LineWidth', 1.5);
%     title(['Joint \alpha_' num2str(i)]);
%     xlabel('Time (s)'); ylabel('Torque (Nm)');
% end
% % 
figure;
for i = 1:n
    subplot(2,1,i);
    plot(tspan, d_mat(:,i), 'k', 'LineWidth', 1.5);
    title(['Joint d' num2str(i)]);
    xlabel('Time (s)'); ylabel('d');
end


function [rho] = performance_poly1(t)
% 保护 beta 在 [0,1] 内
% if t>0
T_p=3;
rho0=2;
p    = 0.3;
s0   = t/T_p;
a=0.02;
s    = min(max(s0^p,0),1);
h    = 1 - 3*s.^2 + 2*s.^3;
dh   = -6*s     + 6*s.^2;                  % dh/ds
dsdt = p * s0.^(p-1) / T_p;                % ds/dt
rho  = a + (rho0 - a)*h;
end

function [rho] = performance_poly2(t)
% 保护 beta 在 [0,1] 内
% if t>0
T_p=3;
rho0=4;
p    = 0.3;
s0   = t/T_p;
a=0.01;
s    = min(max(s0^p,0),1);
h    = 1 - 3*s.^2 + 2*s.^3;
dh   = -6*s     + 6*s.^2;                  % dh/ds
dsdt = p * s0.^(p-1) / T_p;                % ds/dt
rho  = a + (rho0 - a)*h;
end
