function [B, r1, r2, l0_seq, P_m, limb_dir, joint_u_angle_tilt] = deparameterize(p_seq, P_m_nom, limb_dir_nom)
%DEPARAMETERIZE 从标定后的运动学参数 p_seq 反解几何参数（parameterize 的逆）
%   输入：
%     p_seq        — 6×34 运动学参数矩阵（标定后）
%     P_m_nom      — 3×5 名义动平台铰点坐标，提供 x,y 分量（P_m(1:2,:) 不在 p_seq 中）
%     limb_dir_nom — (可选) 5×1 或 5×2 名义支链方向角，用于填充 SPR 支链（第 1 行）
%   输出：
%     B                  — 3×5 基座铰点坐标
%     r1, r2             — 动平台半径（下/上）
%     l0_seq             — 1×5 初始杆长
%     P_m                — 3×5 动平台铰点坐标（z 分量从 p_seq 更新，x,y 沿用 P_m_nom）
%     limb_dir           — 5×2 支链方向角 [θ_base, θ_move]（rad）；
%                           SPR 支链方向角未编码在 p_seq 中，需由 limb_dir_nom 提供
%     joint_u_angle_tilt — U 副倾角（rad），由 4 条 UPS 支链平均得到

    if nargin < 3
        limb_dir_nom = zeros(5, 2);
    end

    P_m = P_m_nom;
    limb_dir = zeros(5, 2);
    limb_dir(1, :) = limb_dir_nom(1, :);  % SPR 支链方向角沿用名义值

    %% SPR 支链（第 1 列组：p_seq(:, 1:6)）
    T1_1 = exp_se3(p_seq(:, 1));
    B(:, 1) = T1_1(1:3, 4);

    T1_5 = exp_se3(p_seq(:, 5));
    l0_seq(1) = T1_5(3, 4);

    T1_p = exp_se3(p_seq(:, 6));
    P_m(3, 1) = T1_p(1, 4);
    r1_from_spr = T1_p(2, 4);

    %% UPS 支链（第 2–5 列组）
    alpha_vals = zeros(4, 1);
    r1_ups = zeros(2, 1);
    r2_ups = zeros(2, 1);

    for i_limb = 2 : 5
        c0 = 7 * (i_limb - 1);  % 该支链 p_seq 起始列（0-based offset）

        % --- B(:, i_limb)：T01 平移分量 ---
        T01 = exp_se3(p_seq(:, c0 + 0));
        B(:, i_limb) = T01(1:3, 4);

        % --- limb_dir(i_limb, 1) 与 joint_u_angle_tilt：T01 旋转矩阵 z 轴 ---
        % T01(:,3) = [sin(α)·cos(θ); sin(α)·sin(θ); cos(α)]
        z01 = T01(1:3, 3);
        alpha_vals(i_limb - 1) = acos(max(min(z01(3), 1), -1));
        limb_dir(i_limb, 1) = atan2(z01(2), z01(1));

        % --- l0_seq(i_limb)：T34 平移分量 ---
        T34 = exp_se3(p_seq(:, c0 + 3));
        l0_seq(i_limb) = T34(3, 4);

        % --- limb_dir(i_limb, 2)：由 T45 的 dif_limb_dir 反推 ---
        % T45(:,3) = [cos(π+diff); sin(π+diff); 0]
        T45 = exp_se3(p_seq(:, c0 + 4));
        z45 = T45(1:3, 3);
        dif = atan2(z45(2), z45(1)) - pi;
        limb_dir(i_limb, 2) = limb_dir(i_limb, 1) - dif;

        % --- P_m(3,i) 与 r：T_p 平移分量 ---
        % SPR/UPS 支链 2,3: r = r1,  t_p = [P_m(3,i); 0; r1]
        % UPS 支链 4,5:       r = r2,  t_p = [P_m(3,i); 0; r2]
        T_p = exp_se3(p_seq(:, c0 + 6));
        P_m(3, i_limb) = T_p(1, 4);
        if i_limb <= 3
            r1_ups(i_limb - 1) = T_p(3, 4);
        else
            r2_ups(i_limb - 3) = T_p(3, 4);
        end
    end

    joint_u_angle_tilt = mean(alpha_vals);
    r1 = mean([r1_from_spr; r1_ups]);
    r2 = mean(r2_ups);

    % 若名义 limb_dir 为单列（仅 θ_base），则输出也裁剪为单列
    if size(limb_dir_nom, 2) == 1
        limb_dir = limb_dir(:, 1);
    end

end
