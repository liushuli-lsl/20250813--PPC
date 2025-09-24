function [w, nu, d_true] = disturbance(t, q, dq, mode, level, band, seed)
% disturbance: 生成并注入扰动信号（仅在 2–3 s 与 5–7 s 作用）
%
% 输入:
%   t      - 当前时间 (s)
%   q,dq   - 系统状态
%   mode   - 'none' | 'matched' | 'unmatched' | 'both'
%   level  - 峰值幅值 (n×1)
%   band   - 基波频率 (rad/s, n×1)
%   seed   - 随机数种子（保证可复现）
%
% 输出:
%   w      - 匹配扰动 (加在加速度通道里, M^{-1}τ 型)
%   nu     - 非匹配扰动 (并入 f(x) 项)
%   d_true - 本次返回的总扰动（已门控）

n = length(q);

% -------- 安全形状处理 --------
if isrow(level), level = level.'; end
if isrow(band),  band  = band.';  end
level = pad_or_clip(level, n);
band  = pad_or_clip(band,  n);

% -------- 只设一次随机种子 --------
persistent seeded
if isempty(seeded) || ~seeded
    rng(seed);
    seeded = true;
end

% -------- 时间门控（2–3 s 与 5–7 s）--------
% gate_type='hard'：矩形门；'smooth'：C^1 平滑门（避免边缘尖峰）
gate_type = 'hard';
win = time_gate(t, gate_type);

% ====== 基本分量（先生成，再乘时间窗） ======
d_sine  = level .* sin(band .* t + (0:n-1)'*pi/4);
% 不再在 t>1 全局加阶跃，改为只在窗口内加一个常值偏置：
d_step  = 0.5 * level;
% 有界噪声
d_noise = 0.2 * level .* (2*rand(n,1) - 1);

% 加权合成 + 时间门控
d_base  = d_sine + d_step + d_noise;
d_true  = win * d_base;       % 只在 2–3 / 5–7 s 生效

% ====== 匹配/非匹配拆分 ======
switch lower(mode)
    case 'matched'
        w  = d_true;                  % 匹配扰动
        nu = zeros(n,1);
    case 'unmatched'
        w  = zeros(n,1);
        % 例：速度相关的非匹配扰动（摩擦/漂移），同样受时间门控
        nu = win * ( 0.1*sign(dq).*level + 0.05*dq.*(level>0) );
    case 'both'
        w  = d_true;
        nu = win * ( 0.1*sign(dq).*level + 0.05*dq.*(level>0) );
    otherwise  % 'none'
        w  = zeros(n,1);
        nu = zeros(n,1);
        d_true = zeros(n,1);
end

end

% ===== 工具：时间门控 =====
function win = time_gate(t, gate_type)
% 返回标量窗函数：在 [2,3] ∪ [5,7] 为 1，其他为 0 或平滑过渡
if strcmpi(gate_type,'smooth')
    % 平滑版本（C^1）：使用软矩形窗和软并集
    w23 = smooth_box(t, 2, 3, 0.05);
    w57 = smooth_box(t, 5, 7, 0.05);
    win = soft_or(w23, w57);
else
    % 硬门
    win = double( (t >= 2 && t <= 3) || (t >= 5 && t <= 7) );
end
end

function y = smooth_box(t, a, b, eps)
% C^1 平滑矩形窗，边缘宽度 eps（秒）
% y≈1 (t∈[a,b])，边缘用 exp(-1/x^2) 平滑
s = @(x) exp(-1./max(x,eps).^2);
left  = s(t - a);
right = s(b - t);
y = (left ./ (left + s(a - t))) .* (right ./ (right + s(t - b)));
end

function z = soft_or(x, y)
% 两个[0,1]窗的“软并集”
z = 1 - (1 - x).*(1 - y);
end

function v = pad_or_clip(v, n)
% 若长度不足，补到 n；若超出，裁剪到 n
m = length(v);
if m < n
    v = [v; repmat(v(end), n-m, 1)];
elseif m > n
    v = v(1:n);
end
end
