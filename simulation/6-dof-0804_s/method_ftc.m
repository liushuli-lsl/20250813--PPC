% 六自由度预定义时间控制器，稳定优化版本 + 动力学模型 + 绘图

clear all;
close all;

%% 控制器参数（减小增益提升稳定性）
Ks = 5;
Kzeta = eye(6);
Ki = 0.2;
D = 0.005 * eye(6);
R = 0.005 * eye(6);
phi = 10;

%% 非线性误差权重
alpha = 0.8; beta = 0.8; p = 0.6; g = 1.1; k = 1.2;
v1 = 0.5; v2 = 0.3; v3 = 0.6; v4 = 1;
rho1 = 1; rho2 = 1; rho3 = 1; rho4 = 1;
d_1 = 0.3 * eye(6);

%% 初始状态
x = [2*pi/4; -2*pi/4;   2*pi/3;  -2*pi/3;   2*pi/4;  -2*pi/4; zeros(6,1)];
zeta = zeros(6,1);
tau = zeros(6,1);
tauH = [60; 60; 40; 30; 20; 15]*2;
tauL =- [60; 60; 40; 30; 20; 15]*2;

%% RBFNN 参数
cij = -5:0.5:5; Node_c = length(cij);
c_c = repmat(cij, 6, 1);
Wc = zeros(Node_c,1);
Sc = zeros(Node_c,1); lrc = 0.05; width_c = 0.2;
c_a = repmat(cij, 12, 1); Node_a = size(c_a,2);
Wa = zeros(Node_a,6);
lra = 1.0; width_a = 2;

%% 仿真设置
Time = 5; step_size = 0.001;
Steps = Time / step_size;
out = zeros(54, floor(Steps/10));
rec = 1;

for count = 2:Steps
    t = (count - 1) * step_size;
    % 轨迹生成
    xd = [0.1*sin(0.5*t)+cos(0.5*t);
           0.1*sin(t)+cos(t);
           0.2*sin(1.5*t)+0.8*cos(t);
           0.3*sin(2*t)+0.7*cos(0.5*t);
           0.1*sin(0.3*t)+0.9*cos(0.2*t);
           0.4*sin(t)+0.6*cos(2*t)];
    dxd = [0.05*cos(0.5*t)-0.5*sin(0.5*t);
           0.1*cos(t)-sin(t);
           0.3*cos(1.5*t)-0.8*sin(t);
           0.6*cos(2*t)-0.35*sin(0.5*t);
           0.03*cos(0.3*t)-0.18*sin(0.2*t);
           0.4*cos(t)-1.2*sin(2*t)];
    ddxd = [-0.025*sin(0.5*t)-0.25*cos(0.5*t);
            -0.1*sin(t)-cos(t);
            -0.45*sin(1.5*t)-0.8*cos(t);
            -1.2*sin(2*t)-0.175*cos(0.5*t);
            -0.009*sin(0.3*t)-0.036*cos(0.2*t);
            -0.4*sin(t)-2.4*cos(2*t)];

    x1 = x(1:6); x2 = x(7:12);
    e1 = x1 - xd;
    e2 = x2 - dxd;

    % Critic 网络
    for i = 1:Node_c
        Sc(i,1) = exp(-(e1 - c_c(:,i))' * (e1 - c_c(:,i)) / width_c^2);
    end
    hatI = Wc' * Sc;
    for i = 1:Node_c
        A(i,1) = -(Sc(i)/phi) + (-2*Sc(i)*(e1 - c_c(:,i))' / width_c^2) * e2;
    end
    dWc = -lrc * (e1'*D*e1 + tau'*R*tau + Wc'*A) * A;

    % 滑模向量 S
    ke1a = (alpha * abs(e1).^(p - 1/(k*v1)) + beta * abs(e1).^(g - 1/(k*v1))).^(k*v1);
    ke1b = k*v1 * (alpha * abs(e1).^(p - 1/(k*v1)) + beta * abs(e1).^(g - 1/(k*v1))).^(k*v1-1) .* ...
             (alpha * (p - 1/(k*v1)) * abs(e1).^(p - 1) + beta * (g - 1/(k*v1)) * abs(e1).^(g - 1));
    Ke1 = diag(ke1a); Ke2 = diag(ke1b);
    S = Ke1 * e1 + tanh(e2) .* abs(e2).^v1;

    % Actor 网络
    za = [e1; e2];
    for i = 1:Node_a
        Sa(i,1) = exp(-(za - c_a(:,i))' * (za - c_a(:,i)) / (width_a)^2);
    end
    Fnn = Wa' * Sa;
    tanh_input = Fnn + Ki * hatI;
    prho = Sa * tanh_input';
    dWa = -lra * prho;

    [M, C, G11] = six_link_dynamics(x1, x2);

    AA1 = -inv(v1)*(Ke1 + Ke2)*tanh(e2) .* abs(e2).^(2-v1);
    BB1 = -inv(v1)*diag(abs(e2).^(1-v1))*(tanh(rho1*S.^v2 + rho2*S.^v3) + zeta + Ks*S);
    tau0 = M*(AA1 + BB1 + ddxd - Fnn - d_1*tanh(S));

    tau = max(min(tau0, tauH), tauL);
    dtau = tau - tau0;
    dzeta = -Kzeta*zeta - tanh(rho3*zeta.^v2 + rho4*zeta.^v3) + S + dtau;

    dx1 = x2;
    dx2 = M \ (tau - C - G11);
    x = x + step_size * [dx1; dx2];
    zeta = zeta + step_size * dzeta;
    Wc = Wc + step_size * dWc;
    Wa = Wa + step_size * dWa;

    if mod(count,10) == 0
        rec = rec + 1;
        out(:,rec) = [x1; x2; xd; dxd; e1; e2; tau; S; Fnn];
    end
end

%% 绘图
T = linspace(0, Time, size(out,2));
figure;
for i = 1:6
    subplot(3,2,i);
    plot(T, out(i,:) - out(12+i,:), 'LineWidth',1.2); hold on;
    title(['Tracking Error e_',num2str(i)]);
    xlabel('Time [s]'); ylabel('Error'); grid on;
end
figure;
for i = 1:6
    subplot(3,2,i);
    plot(T, out(6*6+i,:) , 'LineWidth',1.2); hold on;
    title(['Sliding Surface S_',num2str(i)]);
    xlabel('Time [s]'); ylabel('Value'); grid on;
end
sgtitle('Tracking Error and Sliding Surface Evolution');
