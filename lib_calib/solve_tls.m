function [x, info] = solve_tls(A, b, lambda_ext)
%SOLVE_TLS 总体最小二乘求解 Ax ≈ b（含 Tikhonov 正则化与外部阻尼）
%   参考文献：c02 并联机构的运动学误差建模及参数可辨识性分析_孔令雨 第四章
%
%   组合正则化策略：
%       x = V * diag(s_i / (s_i² - σ_{n+1}² + λ_ext)) * U' * b
%   - λ_ext = 0 且 TLS 适定时：纯 TLS（式 4-18）
%   - λ_ext = σ_{n+1}²：退化为标准最小二乘
%   - λ_ext > σ_{n+1}²：阻尼最小二乘（类似 LM）
%
%   输入：
%     A          — m×n 辨识雅可比 Jp_bar
%     b          — m×1 位姿误差向量
%     lambda_ext — 可选，外部阻尼参数（默认 0，纯 TLS）
%
%   输出：
%     x    — n×1 参数增量
%     info — 结构体 (.method, .eff_rank, .gap, .sigma_np1)

    if nargin < 3
        lambda_ext = 0;
    end

    [Ua, Sa, Va] = svd(A, 'econ');
    sa = diag(Sa);
    n = length(sa);

    % 增广矩阵 SVD → TLS 噪声水平 σ_{n+1}
    [~, Sc, Vc] = svd([A, -b], 'econ');
    sigma_c = diag(Sc);
    sigma_np1 = sigma_c(end);
    gap = sigma_c(end-1) - sigma_c(end);

    % ---- TLS 适定性判定 (Golub & Van Loan) ----
    % TLS 存在唯一解当且仅当 σ_n(A) > σ_{n+1}([A, -b])
    sa_min = sa(end);
    tls_well_posed = (sa_min > sigma_np1 * 1.01);

    % 增广矩阵最后一列右奇异向量 → TLS 原始解
    v_last = Vc(:, end);
    v22 = v_last(end);

    % 组合正则化参数：λ_net = λ_ext - σ_{n+1}²
    % λ_net < 0 → TLS 体制（反阻尼）；λ_net > 0 → 阻尼 LS
    sigma_np1_sq = sigma_np1^2;
    lambda_net = lambda_ext - sigma_np1_sq;

    % 逐分量滤波：d_i = s_i / (s_i² + λ_net)
    % 对 TLS 体制（λ_net < 0），限制最大放大倍数防止数值爆炸
    denom = sa.^2 + lambda_net;
    max_amplify = 100;  % 单分量最多放大 100 倍
    denom = max(denom, sa.^2 / max_amplify);
    d = sa ./ denom;
    x = Va * (d .* (Ua' * b));

    % 分类标记
    if tls_well_posed && lambda_ext < sigma_np1_sq * 0.5 && abs(v22) > 1e-10
        info.method = 'TLS';
    elseif lambda_ext > sigma_np1_sq * 2
        info.method = 'DLS';
    else
        info.method = 'RTLS';
    end
    info.eff_rank = sum(sa.^2 + lambda_net > sa.^2 / max_amplify);
    info.gap = gap;
    info.sigma_np1 = sigma_np1;

end
