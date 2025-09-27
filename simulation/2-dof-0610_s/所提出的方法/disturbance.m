function [w, nu, d_true,tau_d] = disturbance(t, q, dq, mode, level, band, seed)
% disturbance: 在指定时间窗 2–3s 与 5–7s 内施加“阶跃+低频正弦”扰动
%   w  - 匹配扰动（加在加速度/输入通道 g*tau 处）
%   nu - 非匹配扰动（并入 f(x) 漂移）
%   d_true - 实际生成的匹配扰动信号（用于记录/绘图）
%
% 说明：
% - level: 峰值尺度(n×1)，控制阶跃与正弦的幅值基准
% - band : 低频(rad/s, n×1)，建议很小，如 [0.5;0.8]；若为 0 则仅阶跃
% - mode : 'none' | 'matched' | 'unmatched' | 'both'
% - 采用时间门控，仅在 2–3s 与 5–7s 生效

n = length(q);

% ---- 形状安全：统一列向量并裁剪/填充到 n ----
if isrow(level), level = level.'; end
if isrow(band),  band  = band.';  end
level = pad_or_clip(level, n);
band  = pad_or_clip(band,  n);

% ---- 只设一次随机种子 ----
persistent seeded
if isempty(seeded) || ~seeded
    rng(seed);
    seeded = true;
end

% ---- 时间门控（硬门或平滑门二选一）----
gate_type = 'hard';                         % 'hard' | 'smooth'
win = time_gate(t, gate_type)  ;            % 标量 ∈ {0,1} 或 (0,1)

 if t>=11 && t<=13
     level=level*0.15;
%  else
%      level=level;
 end
% ---- 组成：阶跃 + 低频正弦 (+ 少量噪声，可为 0) ----
step_gain = 1.00;                            % 阶跃系数（=level的倍数）
sine_gain = 0.60;                            % 正弦系数
noise_gain = 0.05;                           % 噪声系数（0 关闭）

phase = (0:n-1)' * pi/4;                     % 通道相位错开，便于区分
% 低频正弦：当 band(i)=0 时，此项自动为 0（只剩阶跃）
d_sine  = sine_gain * (level .* sin(band.*t + phase));
d_step  = step_gain * level;                 % 常值偏置（窗口内持续推动）
d_noise = noise_gain * (level .* (2*rand(n,1)-1));

d_base  = d_sine +0;
d_true  = win * d_base;                      % 仅在窗口内有效

% 在窗口内生成力矩阶跃+低频正弦

tau_step = 1.0 * level;                          % 力矩阶跃 (N·m)
tau_sine = 0.6 * (level .* sin(band.*t + phase));% 低频正弦 (N·m)
tau_noise= 0.05* (level .* (2*rand(n,1)-1));     % 少量噪声 (N·m)
tau_d_base = tau_step + tau_sine + tau_noise;
tau_d = win * tau_d_base;
% ---- 匹配/非匹配通道拆分 ----
switch lower(mode)
  case 'matched'
    w  = d_true;
    nu = zeros(n,1);
  case 'unmatched'
    w  = zeros(n,1);
    % 速度相关非匹配项（科氏/摩擦等效）：同样受门控
    nu = win * ( 0.10*sign(dq).*level + 0.05*dq.*(level>0) );
  case 'both'
         w  = d_true;
       tau_d=tau_d;
       nu = win * ( 0.10*sign(dq).*level + 0.05*dq );
  otherwise % 'none'
    w  = zeros(n,1);
    nu = zeros(n,1);
    d_true = zeros(n,1);
end
end

% ===== 辅助函数 =====
function v = pad_or_clip(v, n)
m = length(v);
if m < n
    v = [v; repmat(v(end), n-m, 1)];
elseif m > n
    v = v(1:n);
end
end

function win = time_gate(t, gate_type)
if strcmpi(gate_type,'smooth')
    % 平滑门（C^1）：2–3s 与 5–7s，边缘 0.05s 过渡
    w23 = smooth_box(t, 2, 3, 0.05);
    w57 = smooth_box(t, 5, 7, 0.05);
    win = 1 - (1 - w23).*(1 - w57);  % 软并集
else
    % 硬门
    win = double( (t>=4 && t<=5) || (t>=11 && t<=13) );
end
end

function y = smooth_box(t, a, b, eps)
% C^1 平滑矩形窗，边缘宽 eps（秒）
s = @(x) exp(-1./max(x,eps).^2);
left  = s(t - a);  right = s(b - t);
y = (left ./ (left + s(a - t))) .* (right ./ (right + s(t - b)));
end




