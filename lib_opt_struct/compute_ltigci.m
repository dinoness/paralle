function [LTI, GCI, OTI, ITI] = compute_ltigci(T_ref, B, l0_seq, P_m, p_seq)
% 根据 force_analysis_screw.m 计算指定姿态下的 LTI、GCI、OTI 和 ITI
%
% 输入：
%   T_ref  : 4×4 位姿变换矩阵
%   B      : 3×5 基座关节坐标
%   l0_seq : 5×1 初始支链长度
%   P_m    : 3×5 动平台关节在平台坐标系下的坐标
%   p_seq  : 6×34 旋量参数（由 parameterize 预计算）
%
% 输出：
%   LTI    : 局部传递指标
%   GCI    : 全局约束指标
%   OTI    : 输出传递指标
%   ITI    : 输入传递指标

    zeta_r = [0;0;1;0;0;0];
    zeta_p = [0;0;0;0;0;1];
    Omega = [zeros(3,3) eye(3); eye(3) zeros(3,3)];

    joint_q0 = keni_sol_inverse(T_ref, B, l0_seq, P_m, p_seq);
    P = T_ref(1:3, 4) + T_ref(1:3, 1:3)*P_m;

    U = zeros(6, 6, 5);
    ST = zeros(6, 5);
    SC = zeros(6, 1);
    SI = zeros(6, 5);
    SO = zeros(6, 5);
    Delta_SO = zeros(6, 1);

    % SPR 支链
    U1 = zeros(6, 6);
    T1 = zeros(4,4,6);
    for i_joint = 1 : 6
        T1(:,:,i_joint) = exp_se3(p_seq(:, i_joint));
    end

    T_temp = eye(4);
    T_zeta = eye(4);
    for i_joint = 1 : 5
        T_temp = T_temp * T_zeta * T1(:,:,i_joint);

        if i_joint ~= 4
            T_zeta = exp_se3(zeta_r * joint_q0(i_joint,1));
        else
            T_zeta = exp_se3(zeta_p * joint_q0(i_joint,1));
        end

        if i_joint ~= 4
            U1(:, i_joint) = adjoint_m(T_temp) * zeta_r;
        else
            U1(:, i_joint) = adjoint_m(T_temp) * zeta_p;
        end
    end
    U(:, :, 1) = U1;
    SI(:, 1) = U1(:, 4);
    SC(:, 1) = null(U1(:, 1:5)' * Omega);
    WA1 = Omega * U1(:, 1:5) / (U1(:, 1:5)' * U1(:, 1:5));
    TR1 = null(WA1' * Omega);
    ST(:, 1) = [U1(4:6,4); cross(B(:, 1), U1(4:6,4))];

    % UPS 支链
    T2 = zeros(4,4,7);
    for i_limb = 2 : 5
        U2 = zeros(6, 6);
        for i_joint = 1 : 7
            T2(:,:,i_joint) = exp_se3(p_seq(:, (7*(i_limb-1) - 1) + i_joint));
        end

        T_temp = eye(4);
        T_zeta = eye(4);
        for i_joint = 1 : 6
            T_temp = T_temp * T_zeta * T2(:,:,i_joint);

            if i_joint ~= 3
                T_zeta = exp_se3(zeta_r * joint_q0(i_joint,i_limb));
            else
                T_zeta = exp_se3(zeta_p * joint_q0(i_joint,i_limb));
            end

            if i_joint ~= 3
                U2(:, i_joint) = adjoint_m(T_temp) * zeta_r;
            else
                U2(:, i_joint) = adjoint_m(T_temp) * zeta_p;
            end
        end

        U(:, :, i_limb) = U2;
        ST(:, i_limb) = [U2(4:6,3); cross(B(:,i_limb), U2(4:6,3))];
        SI(:, i_limb) = U2(:, 3);
    end

    % 计算输出运动旋量 SO，与受限运动旋量 Delta_SO
    for i = 1 : 5
        idx = [1:i-1, i+1:5];
        W = [SC, ST(:, idx)];
        SO(:, i) = null(W' * Omega);
        if ST(:, i)' * Omega * SO(:, i) < 0
            SO(:, i) = -SO(:, i);
        end
    end

    Delta_SO(:, 1) = null(ST' * Omega);

    % ITI / OTI
    lambda = zeros(5, 1);
    eta = zeros(5, 1);
    for i = 1: 5
        lambda(i) = screw_efficiency(ST(:, i), SI(:, i), 'p_cstr', P(:, i));
        eta(i) = screw_efficiency(ST(:, i), SO(:, i), 'p_cstr', P(:, i));
    end

    ITI = min(lambda);
    OTI = min(eta);
    LTI = min(ITI, OTI);

    % ICI / OCI
    zeta = screw_efficiency(SC(:, 1), TR1, 'p_cstr', B(:, 1));
    kappa = screw_efficiency(SC(:, 1), Delta_SO(:, 1), 'p_cstr', B(:, 1));
    GCI = min(zeta, kappa);

end
