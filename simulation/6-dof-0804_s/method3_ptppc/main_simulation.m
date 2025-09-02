clearvars;          % 清除工作区变量
clear functions;    % 清除所有函数的 persistent 缓存
close all;          % 关闭所有图窗
clc;                % 清空命令窗口


%% 1) 仿真设置
n=6;  %自由度
T_total = 20; 
N       = 20000;            % 点数
tspan   = linspace(0, T_total, N);
dt      = tspan(2) - tspan(1);

%% 2) 初始状态
% q00 = [
%           pi,     pi;         % 偏离参考较大
%     pi/2,     pi/2;        % 中等伸展                     T=3.5,kI=[12, 6]*10;
% 0,    pi/2;        % 十字形对称
%  pi/4,     pi/4;        % 十字形对称
%       -pi/2,    -pi/2;        % 十字形对称
%        -3*pi/4,     -3*pi/4;      % 偏离参考较大
%           0,        0;          % 折叠（近奇异）
% ];
q00 = [

    pi/4,    -pi/8,     pi/6,     -pi/4,     pi/8,     -pi/8;    %可以运行
    1.12, 1.08, 0.95, 0.8, 0.95, 0.67;  %k可以运行
         0,  -0.52,   0,  1.79,  0,  0.25;   % Case 1: 30°, -30°, 45°, -45°, 60°, -60°
    2*pi/4,  -2*pi/4,   2*pi/3,  -2*pi/3,   2*pi/4,  -2*pi/4;   % Case 2: 90°, -90°,120°,-120°,135°,-135°
    2.62,  -2.62,   2.36,  -2.36,   2.62,  -2.62;    % Case 3:150°,-150°,135°,-135°,150°,-150°
    pi, 1.0, 0.8, 0.7, 0.9, 0.6;
];
q0 = q00(4, :)';   % q0 = [0; 0]  2,3,4,5
dq0 = [0;0;0;0;0;0];                  %初始誤差不能爲0，否則會產生歧義
% q0  = zeros(6,1);
% dq0 = zeros(6,1);
zeta1   = zeros(n,1);
zeta2   = zeros(n,1);
e1_int  = zeros(n,1);
alpha_bar = zeros(n,1);
d1 = zeros(n,1);
d2 = zeros(n,1);
c=[2;2];
x0  = [q0; dq0;zeta1; zeta2;e1_int;alpha_bar;c;d1;d2];



%% 3) 预分配存储

x       = zeros(8*n+2, N);
q_use   = zeros(N,n);
dq_use  = zeros(N,n);
qd_mat  = zeros(N,n);
dqd_mat = zeros(N,n);
tau_mat = zeros(N,n);
alpha_mat = zeros(N,n);
W_hat_mat = zeros(N,n);
zeta1_mat= zeros(N,n);
zeta2_mat= zeros(N,n);

x(:,1) = x0;
q_use(1,:)  = q0.';
dq_use(1,:) = dq0.';

%% 4) 固定步长 RK4 主循环
for k = 1:N-1
    t = tspan(k);
    % 解算一次 RK4
  [k1, tau1, alpha1,W_hat1] = controller_ptc(t,           x(:,k));
[k2, tau2, alpha2,W_hat2] = controller_ptc(t+dt/2,      x(:,k)+dt/2*k1);
[k3, tau3, alpha3,W_hat3] = controller_ptc(t+dt/2,      x(:,k)+dt/2*k2);
[k4, tau4, alpha4,W_hat4] = controller_ptc(t+dt,        x(:,k)+dt*k3);

    x(:,k+1) = x(:,k) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;

    % 存储关节、参考、控制
    q   = x(1:1*n,k+1);
    dq  = x(n+1:2*n,k+1);
    zeta1= x(2*n+1:3*n,k+1);
     zeta2= x(3*n+1:4*n,k+1);

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

    q_use(k+1,:)   = q.';
    dq_use(k+1,:)  = dq.';
    qd_mat(k,:)  = qd.';
    dqd_mat(k,:) = dqd.';

%     tau = controller_ptc(t+dt, q, dq);
    tau_mat(k+1,:) = tau1.';
     alpha_mat(k+1,:) = alpha1.';
     zeta1_mat(k+1,:) = zeta1.';
     zeta2_mat(k+1,:) = zeta2.';
%      W_hat_mat(k+1,:) = W_hat1.';
end

e_q  = q_use   - qd_mat;
e_dq = dq_use  - dqd_mat;

% === 在 main_simulation.m 的绘图部分后面添加 ===
% alpha_use  = ALPHA_BAR; 


[rho1] = arrayfun(@(tt) performance_poly1(tt),tspan);
[rho2] = arrayfun(@(tt) performance_poly2(tt),tspan);

save('method3.mat', 'tspan', 'e_q', 'e_dq', 'tau_mat', 'qd_mat', 'q_use','dqd_mat', 'dq_use','rho1','rho2',"zeta1_mat",'zeta2_mat',"alpha_mat");
figure;
for i = 1:n
    subplot(3,2,i);
    plot(tspan, e_q(:,i), 'b', 'LineWidth', 1.5); hold on;
    plot(tspan, rho1, 'k--', 'LineWidth',1.5);hold on;
      plot(tspan, -rho1, 'k--', 'LineWidth',1.5);
    yline(0, 'k--');
    xline(3, 'r--', 'LineWidth', 1.2);
    title(['Tracking Error e_' num2str(i)]);
    xlabel('Time (s)'); ylabel('e_i (rad)');
%     legend('Error', 'Zero Line', 'T_p');
end


% Plot error curves with T_p marker
figure;
for i = 1:n
    subplot(3,2,i);
    plot(tspan, e_dq(:,i), 'b', 'LineWidth', 1.5); hold on;
      plot(tspan, rho2, 'k--', 'LineWidth',1.5);hold on;
      plot(tspan, -rho2, 'k--', 'LineWidth',1.5);
    yline(0, 'k--');
    xline(3, 'r--', 'LineWidth', 1.2);
    title(['Tracking Error e_' num2str(i)]);
    xlabel('Time (s)'); ylabel('e_i (rad)');
%     legend('Error', 'Zero Line', 'T_p');
end





% 
% 
% % Plot joint position vs reference
figure;
for i = 1:n
    subplot(3,2,i);
    plot(tspan, q_use(:,i), 'b', 'LineWidth', 1.5); hold on;
    plot(tspan, qd_mat(:,i), 'r--', 'LineWidth', 1.2);
    title(['Joint q_' num2str(i)]);
    legend('Actual', 'Reference');
    xlabel('Time (s)'); ylabel('Position (rad)');
end




% Plot joint velocity vs reference
figure;
for i = 1:n
    subplot(3,2,i);
    plot(tspan, dq_use(:,i), 'b', 'LineWidth', 1.5); hold on;
    plot(tspan, dqd_mat(:,i), 'r--', 'LineWidth', 1.2);
    title(['Joint dq_' num2str(i)]);
    legend('Actual', 'Reference');
    xlabel('Time (s)'); ylabel('Velocity (rad/s)');
end
% % 
% % Plot control torques
figure;
for i = 1:n
    subplot(3,2,i);
    plot(tspan, tau_mat(:,i), 'k', 'LineWidth', 1.5);
    title(['Joint \tau_' num2str(i)]);
    xlabel('Time (s)'); ylabel('Torque (Nm)');
end
figure;
for i = 1:n
    subplot(3,2,i);
    plot(tspan, alpha_mat(:,i), 'k', 'LineWidth', 1.5);
    title(['Joint \alpha_' num2str(i)]);
    xlabel('Time (s)'); ylabel('Torque (Nm)');
end
% 
figure;
for i = 1:n
    subplot(3,2,i);
    plot(tspan, zeta1_mat(:,i), 'k', 'LineWidth', 1.5);
    title(['Joint \zeta_1' num2str(i)]);
    xlabel('Time (s)'); ylabel('Torque (Nm)');
end
figure;
for i = 1:n
    subplot(3,2,i);
    plot(tspan, zeta2_mat(:,i), 'k', 'LineWidth', 1.5);
    title(['Joint \zeta_2' num2str(i)]);
    xlabel('Time (s)'); ylabel('Torque (Nm)');
end


function [rho] = performance_poly1(t)
% 保护 beta 在 [0,1] 内
% if t>0
T_p=2.5;
rho0=2*pi/3;
p    = 0.6;
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
T_p=2.5;
rho0=pi;
p    = 0.6;
s0   = t/T_p;
a=0.01;
s    = min(max(s0^p,0),1);
h    = 1 - 3*s.^2 + 2*s.^3;
dh   = -6*s     + 6*s.^2;                  % dh/ds
dsdt = p * s0.^(p-1) / T_p;                % ds/dt
rho  = a + (rho0 - a)*h;
end
