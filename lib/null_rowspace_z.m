function [N, R, Q] = null_rowspace_z(xi)
% 根据论文第58-59页方法，计算 Z_xi 的零空间和行空间标准正交基
% 输入：
%   xi - 6x1 旋量，前3个为角速度 omega，后3个为线速度 v
% 输出：
%   N - 6 x dim(null) 矩阵，列向量是零空间的标准正交基
%   R - 6 x rank(Z)   矩阵，列向量是行空间的标准正交基
%   Q - 6 x 6 正交矩阵，Q = [N, R]

    omega = xi(1:3);
    v = xi(4:6);

    % 判断关节类型
    is_revolute = (norm(omega) > 1e-10);   % 转动关节
    is_prismatic = (~is_revolute) && (norm(v) > 1e-10); % 移动关节

    if is_revolute
        % 转动关节：零空间维数 = 2
        % 基向量1: xi 本身
        xi_vec = xi;
        % 基向量2: [0; omega]   (论文中的 ξ_perp)
        xi_perp = [0; 0; 0; omega(1); omega(2); omega(3)];
        N_init = [xi_vec, xi_perp];
        % 对初始基进行 Gram-Schmidt 正交化，得到零空间标准正交基
        N = gram_schmidt(N_init);
        % 行空间基：求 N 的左零空间（即 Z_xi 的行空间）
        % 注意：Z_xi 的行空间与 N 的左零空间正交
        Z = ad_matrix(xi);   % 构造 Z_xi
        % 取 Z 的行空间，可用 null(N') 或 orth(Z')
        R_init = orth(Z');
        % 对行空间基进行 Gram-Schmidt 正交化（可选）
        R = gram_schmidt(R_init);
    elseif is_prismatic
        % 移动关节：零空间维数 = 4
        % 基向量1-3: [0; e_i]  (i=1,2,3)
        e1 = [0;0;0;1;0;0];
        e2 = [0;0;0;0;1;0];
        e3 = [0;0;0;0;0;1];
        % 基向量4: [v; 0]
        % 实际上移动关节的旋量坐标 ξ = [0; v]，本身也在零空间中，但论文给出的四个基为：
        % [0;e1], [0;e2], [0;e3], [v;0] （见式(3-8)）
        xi_vec = [v; 0;0;0];   % [0; v]
        N_init = [e1, e2, e3, xi_vec];
        % Gram-Schmidt 正交化
        N = gram_schmidt(N_init);
        % 行空间基：求 Z_xi 的行空间，维数 = 2
        Z = ad_matrix(xi);
        R_init = orth(Z');
        R = gram_schmidt(R_init);
    else
        % 一般旋量（理论上不会出现，但作为通用回退）
        Z = ad_matrix(xi);
        N = null(Z);           % MATLAB 自动正交化
        R = orth(Z');          % 行空间标准正交基
    end

    % 组合成正交矩阵 Q
    Q = [N, R];
end

function Z = ad_matrix(xi)
% 构造伴随矩阵 Z_xi
    omega = xi(1:3);
    v = xi(4:6);
    omega_hat = skew(omega);
    v_hat = skew(v);
    Z = [omega_hat, zeros(3);
         v_hat,     omega_hat];
end

function S = skew(v)
% 构造反对称矩阵
    S = [ 0,   -v(3),  v(2);
          v(3), 0,    -v(1);
         -v(2), v(1),  0];
end

function Q = gram_schmidt(V)
% Gram-Schmidt 正交化（列向量）
% 输入 V：n x k 矩阵，列向量线性无关
% 输出 Q：n x k 矩阵，列向量标准正交
    [n, k] = size(V);
    Q = zeros(n, k);
    for i = 1:k
        q = V(:, i);
        for j = 1:i-1
            q = q - (Q(:, j)' * q) * Q(:, j);
        end
        norm_q = norm(q);
        if norm_q > 1e-12
            Q(:, i) = q / norm_q;
        else
            error('列向量线性相关，无法正交化');
        end
    end
end