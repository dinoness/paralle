function [U, V_prep, V] = calib_iter_row_space_matrix(xi_seq)
% inupt ：关节零位全局坐标 1-6,7-13,14-20,21-27,28-34
% output：行空间矩阵

U = cell(1, 5);
V = cell(1, 5);
V_prep = cell(1, 5);

% SPR支链
U1_blk = cell(1, 6);
Xi_bar1_blk = cell(1, 5);

for i_joint = 1 : 5
    Xi_bar1_blk{i_joint} = xi_seq(:, i_joint);
    [~, R] = null_rowspace_z(xi_seq(:, i_joint));
    U1_blk{i_joint} = R;
end

U1_blk{6} = eye(6);
U{1} = blkdiag(U1_blk{:});

Xi_bar1_blk(4) = [];  % 删除移动副
Xi_bar1 = blkdiag(Xi_bar1_blk{:});
Xi_bar1 = [Xi_bar1(1:3*6, :);
           zeros(6, size(Xi_bar1, 2));
           Xi_bar1(3*6+1:end, :);
           zeros(6, size(Xi_bar1, 2))];
Delta_5 = construct_delta(5);  % 实际是36*36，5是因为只有5个关节
V{1} = U{1}' * Delta_5 * Xi_bar1;
V_prep{1} = null(V{1}');


% UPS支链

for i_limb = 2 : 5
    U2_blk = cell(1, 7);
    Xi_bar2_blk = cell(1, 6);
    for i_joint = 1 : 6
        Xi_bar2_blk{i_joint} = xi_seq(:, 7*(i_limb-1) + (i_joint-1));
        [~, R] = null_rowspace_z(xi_seq(:, 7*(i_limb-1) + (i_joint-1)));
        U2_blk{i_joint} = R;
    end

    U2_blk{7} = eye(6);
    U{i_limb} = blkdiag(U2_blk{:});

    Xi_bar2_blk(3) = [];  % 删除移动副
    Xi_bar2 = blkdiag(Xi_bar2_blk{:});
    Xi_bar2 = [Xi_bar2(1:2*6, :);
                zeros(6, size(Xi_bar2, 2));
                Xi_bar2(2*6+1:end, :);
                zeros(6, size(Xi_bar2, 2))];
    
    
    Delta_6 = construct_delta(6);  % 实际是42*42，6是因为只有6个关节
    V{i_limb} = U{i_limb}' * Delta_6 * Xi_bar2;
    V_prep{i_limb} = null(V{i_limb}');
end


end

% 怎么验证U是对的





function Delta = construct_delta(n_i)
    % CONSTRUCT_DELTA 构造分块下三角单位矩阵
    %   输入：
    %       n_i - 非负整数，决定矩阵的维数
    %   输出：
    %       Delta - 大小为 6*(n_i+1) 的分块下三角矩阵

    % 参数检查
    if n_i < 0 || floor(n_i) ~= n_i
        error('n_i 必须为非负整数');
    end

    % 子块数量
    N = n_i + 1;
    
    % 方法：Kronecker 乘积
    % tril(ones(N)) 生成 N×N 下三角全1矩阵（主对角及以下为1，以上为0）
    % 用 eye(6) 替换每个元素1，用 0 替换每个元素0，即得所需分块矩阵
    Delta = kron(tril(ones(N)), eye(6));
end

