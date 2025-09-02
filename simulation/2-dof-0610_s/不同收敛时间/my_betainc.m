function I = my_betainc(x, a, b)
% my_betainc: 计算正则化不完全Beta函数 betainc(x,a,b) 的数值近似
% 用于替代 Statistics Toolbox 中的 betainc 函数

% 保证 x 为行向量
x = x(:)';
I = zeros(size(x));

for i = 1:length(x)
    xi = x(i);
    if xi <= 0
        I(i) = 0;
    elseif xi >= 1
        I(i) = 1;
    else
        % 定义被积函数 f(t) = t^(a-1) * (1 - t)^(b-1)
        f = @(t) t.^(a-1) .* (1 - t).^(b-1);
        % 使用 quadgk 或 integral 计算不完全Beta积分
        numerator = integral(f, 0, xi, 'RelTol',1e-8,'AbsTol',1e-12);
        denominator = beta(a, b);  % 完全Beta函数
        I(i) = numerator / denominator;
    end
end
end
