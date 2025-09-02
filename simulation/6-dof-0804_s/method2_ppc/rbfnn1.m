function [phi,W] = rbfnn(x_in, e2, rho)
    % RBFNN estimate and update for 6-DOF robot
    % x_in: [q; dq] state (12x1)
    % e2: velocity error (6x1)
    % rho: performance bound scalar
n=6;
    persistent C_rbfnn D_rbfnn W_hat

    N = 30;
    state_dim = 2*n;
    if isempty(C_rbfnn)
        C_rbfnn = rand(state_dim, N) * 2 - 1;
        D_rbfnn = 2*ones(N,1);
       W_hat = cell(6,1);
for i = 1:6
    W_hat{i} = zeros(N, 1);
end
    end

    varrho = 10; sigma = 3; step = 0.001;
    phi = zeros(6,1);
    for i = 1:6
        Q_i = zeros(N,1);
        for j = 1:N
            Q_i(j) = exp(-norm(x_in - C_rbfnn(:,j))^2) /D_rbfnn(j)^2;
        end
        phi(i) = W_hat{i}' * Q_i;
        W_hat{i} = W_hat{i} + step *(varrho * Q_i * (e2(i) / (rho^2 - e2(i)^2)) - sigma * W_hat{i});
        W=W_hat;
    end
end
