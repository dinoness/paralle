function delta_p_seq = restore_full_param_increment(pk, U, V_prep, Lambda, Ap)
% 根据缩减后的独立参数增量 pk，恢复原始运动学参数增量 delta_p_seq
% 恢复关系与式(3-28)一致：
%   delta_p_i = A_{p,i}^{-1} * Lambda_i^{-1} * Delta_{n_i}^{-1} * U_{xi,i} * Vbar_perp_i * delta_pbar_i
% 其中 Vbar_perp_i 已在 calib_iter_row_space_matrix 中通过 null(V_i') 构造，
% 因而冗余参数分量已被显式消去，不需要再额外拼接 [delta_pbar_i; 0] 后乘 M_i。

    n_limb = numel(U);
    delta_p_cells = cell(1, n_limb);
    idx_start = 1;
    total_cols = 0;

    for i_limb = 1 : n_limb
        reduced_dim_i = size(V_prep{i_limb}, 2);
        idx_end = idx_start + reduced_dim_i - 1;
        if idx_end > numel(pk)
            error('pk 的长度与各支链缩减参数维数不一致');
        end

        delta_pbar_i = pk(idx_start : idx_end);
        delta_rho_i = U{i_limb} * (V_prep{i_limb} * delta_pbar_i);

        num_blocks_i = size(Ap{i_limb}, 1) / 6;
        Delta_i = construct_delta_matrix(num_blocks_i);

        delta_p_i = Ap{i_limb} \ (Lambda{i_limb} \ (Delta_i \ delta_rho_i));
        delta_p_cells{i_limb} = reshape(delta_p_i, 6, []);

        total_cols = total_cols + size(delta_p_cells{i_limb}, 2);
        idx_start = idx_end + 1;
    end

    if idx_start ~= numel(pk) + 1
        error('pk 中仍有未分配的参数分量');
    end

    delta_p_seq = zeros(6, total_cols);
    col_start = 1;
    for i_limb = 1 : n_limb
        col_end = col_start + size(delta_p_cells{i_limb}, 2) - 1;
        delta_p_seq(:, col_start : col_end) = delta_p_cells{i_limb};
        col_start = col_end + 1;
    end
end


function Delta = construct_delta_matrix(num_blocks)
% 构造 Δ 矩阵，其中 num_blocks = n_i + 1
% 输出矩阵大小为 6*num_blocks × 6*num_blocks

    if num_blocks <= 0 || floor(num_blocks) ~= num_blocks
        error('num_blocks 必须为正整数');
    end

    Delta = kron(tril(ones(num_blocks)), eye(6));
end