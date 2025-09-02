clearvars;          % 清除工作区变量
clear functions;    % 清除所有函数的 persistent 缓存
close all;          % 关闭所有图窗
clc;                % 清空命令窗口


%% 1) 仿真设置
n=2;  %自由度
T_total = 10; 
N       = 10000;            % 点数
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
q0 = q00(2, :)';   % q0 = [0; 0]  2,3,4,5
dq0 = [0;0];                  %初始誤差不能爲0，否則會產生歧義
% q0  = zeros(6,1);
% dq0 = zeros(6,1);
zeta1   = zeros(n,1);
zeta2   = zeros(n,1);
e1_int  = zeros(n,1);
alpha_bar = zeros(n,1);
d1 = zeros(n,1);
d2 = zeros(n,1);
x0  = [q0; dq0;zeta1; zeta2;e1_int;alpha_bar;d1;d2];



%% 3) 预分配存储

x       = zeros(8*n, N);
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
    q   = x(1:2,k+1);
    dq  = x(3:4,k+1);
    zeta1= x(5:6,k+1);
     zeta2= x(7:8,k+1);

qd = [0.1*sin(0.5*t) + cos(0.5*t);0.1*sin(t) + cos(t)];
dqd = [0.05*cos(0.5*t)-0.5*sin(0.5*t); 0.1*cos(t)-sin(t)];

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

save('x2.mat', 'tspan', 'e_q', 'e_dq', 'tau_mat', 'qd_mat', 'q_use','dqd_mat', 'dq_use','rho1','rho2',"zeta1_mat",'zeta2_mat',"alpha_mat");
figure;
for i = 1:n
    subplot(2,1,i);
    plot(tspan, e_q(:,i), 'b', 'LineWidth', 1.5); hold on;
    plot(tspan, rho1, 'k--', 'LineWidth',1.5);hold on;
      plot(tspan, -rho1, 'k--', 'LineWidth',1.5);
    yline(0, 'k--');
    xline(3, 'r--', 'LineWidth', 1.2);
    title(['Tracking Error e_' num2str(i)]);
    xlabel('Time (s)'); ylabel('e_i (rad)');
    legend('Error', 'Zero Line', 'T_p');
end


% Plot error curves with T_p marker
figure;
for i = 1:n
    subplot(2,1,i);
    plot(tspan, e_dq(:,i), 'b', 'LineWidth', 1.5); hold on;
      plot(tspan, rho2, 'k--', 'LineWidth',1.5);hold on;
      plot(tspan, -rho2, 'k--', 'LineWidth',1.5);
    yline(0, 'k--');
    xline(3, 'r--', 'LineWidth', 1.2);
    title(['Tracking Error e_' num2str(i)]);
    xlabel('Time (s)'); ylabel('e_i (rad)');
    legend('Error', 'Zero Line', 'T_p');
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
figure;
for i = 1:n
    subplot(2,1,i);
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
    subplot(2,1,i);
    plot(tspan, tau_mat(:,i), 'k', 'LineWidth', 1.5);
    title(['Joint \tau_' num2str(i)]);
    xlabel('Time (s)'); ylabel('Torque (Nm)');
end
figure;
for i = 1:n
    subplot(2,1,i);
    plot(tspan, alpha_mat(:,i), 'k', 'LineWidth', 1.5);
    title(['Joint \alpha_' num2str(i)]);
    xlabel('Time (s)'); ylabel('Torque (Nm)');
end
% 
figure;
for i = 1:n
    subplot(2,1,i);
    plot(tspan, zeta1_mat(:,i), 'k', 'LineWidth', 1.5);
    title(['Joint \zeta_1' num2str(i)]);
    xlabel('Time (s)'); ylabel('Torque (Nm)');
end
figure;
for i = 1:n
    subplot(2,1,i);
    plot(tspan, zeta2_mat(:,i), 'k', 'LineWidth', 1.5);
    title(['Joint \zeta_2' num2str(i)]);
    xlabel('Time (s)'); ylabel('Torque (Nm)');
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
