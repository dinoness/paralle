function [U, active_idx, passive_idx] = spr4ups_build_limb_twists(p_seq, joint_q, i_limb, zeta_r, zeta_p)
%SPR4UPS_BUILD_LIMB_TWISTS 计算第 i_limb 条支链在当前位形下的关节运动旋量
% 输出：
%   U           : 6xn_joints，每列对应一个关节运动旋量 [omega; v]
%   active_idx  : 主动P副在 U 中的列号
%   passive_idx : 被动关节列号

    if i_limb == 1
        % SPR：5个实际关节，6个定长转移块
        trans_idx = 1 : 6;
        joint_type = {'R', 'R', 'R', 'P', 'R'};
        active_idx = 4;
    else
        % UPS：6个实际关节，7个定长转移块
        trans_idx = (7 * (i_limb - 1)) : (7 * (i_limb - 1) + 6);
        joint_type = {'R', 'R', 'P', 'R', 'R', 'R'};
        active_idx = 3;
    end

    n_joint = numel(joint_type);
    n_trans = numel(trans_idx);

    T_const = zeros(4, 4, n_trans);
    for k = 1 : n_trans
        T_const(:, :, k) = exp_se3(p_seq(:, trans_idx(k)));
    end

    U = zeros(6, n_joint);
    T_prefix = eye(4);
    T_zeta = eye(4);

    for j = 1 : n_joint
        % 到当前关节局部坐标系的变换
        T_prefix = T_prefix * T_zeta * T_const(:, :, j);

        if joint_type{j} == 'R'
            zeta = zeta_r;
        else
            zeta = zeta_p;
        end

        U(:, j) = adjoint_m(T_prefix) * zeta;
        T_zeta = exp_se3(zeta * joint_q(j, i_limb));
    end

    passive_idx = setdiff(1:n_joint, active_idx);
end
