function [ M, C, G] = two_link_dynamics(q,dq)
% TWO_LINK_DYNAMICS  二自由度机械臂刚体动力学
%   [dx, M, C, G] = two_link_dynamics(t, x, tau)
%   输入:
%       t   - 当前时间（可选，用于时变参数或轨迹）
%       x   - 状态向量 [q1; q2; dq1; dq2]
%       tau - 输入扭矩 [tau1; tau2]
%   输出:
%       dx  - 状态导数 [dq1; dq2; ddq1; ddq2]
%       M   - 质量矩阵 2×2
%       C   - Coriolis/Centrifugal 矩阵 2×2
%       G   - 重力向量 2×1

    %%—— 系统参数 ——%%
    m1   = 2;       m2   = 0.85;
    L1   = 0.35;    L2   = 0.31;
    Lc1  = 0.5*L1;  Lc2  = 0.5*L2;   % 质心到关节距离
    I1   = 0.25*m1*L1^2;
    I2   = 0.25*m2*L2^2;
    g    = 9.8;

    %%—— 状态分量 ——%%
    q1  = q(1);
    q2  = q(2);
    dq1 = dq(1);
    dq2 = dq(2);

    %%—— 质量矩阵 M ——%%
    r1 = m1*Lc1^2 + m2*L1^2 + I1;
    r2 = m2*Lc2^2 + I2;
    r3 = m2*L1*Lc2;
    M11 = r1 + r2 + 2*r3*cos(q2);
    M12 = r2 + r3*cos(q2);
    M21 = M12;
    M22 = r2;
    M   = [M11, M12; M21, M22];

    %%—— Coriolis/Centrifugal 矩阵 C ——%%
    h = -r3*sin(q2);
    C = [ h*dq2,  h*(dq1 + dq2);
         -h*dq1,        0       ];

    %%—— 重力向量 G ——%%
    G1 = (m1*Lc1 + m2*L1)*g*cos(q1) + m2*Lc2*g*cos(q1 + q2);
    G2 = m2*Lc2*g*cos(q1 + q2);
    G  = [G1; G2];

    %%—— 状态导数 ——%%
    % 方程: M*ddq + C*dq + G = tau
%     ddq = M \ (tau - C*[dq1;dq2] - G);
% 
%     dx = [dq1; dq2; ddq];

end
