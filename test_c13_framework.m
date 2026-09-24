% c13 标定框架验证：
%   1) dq_fks_spr4ups（对偶四元数正解）与 keni_sol_forward（POE 正解）一致性
%   2) build_emm_j1_spr4ups / build_emm_j3 解析 EMM 与有限差分一致性
clear
path_add();

basic_paras = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'm');
unit_para = basic_paras.unit_para;
kin.B = basic_paras.B;
kin.P_m = basic_paras.P_m;
kin.l0 = basic_paras.l0_seq;
kin.d = 0;   % 非理想约束偏置（第 1 节与 POE 正解对照须取 0）
[p_seq, ~] = parameterize(basic_paras.limb_dir, basic_paras.B, basic_paras.r1, ...
    basic_paras.r2, basic_paras.l0_seq, basic_paras.P_m, basic_paras.joint_u_angle_tilt);

% 任取两个理论位姿
pos_csv = readmatrix(fullfile('calib_p', 'pos_seq.csv'));
Pos_ref_seq = pos_csv.';
Pos_ref_seq(1:3, :) = Pos_ref_seq(1:3, :) * unit_para;
Pos_ref_seq(4:5, :) = deg2rad(Pos_ref_seq(4:5, :));

%% 1) 正解一致性
fprintf('===== 1) DQ-FKS vs keni_sol_forward =====\n');
for im = [1, round(size(Pos_ref_seq,2)/2)]
    T_nom = pos2trans(Pos_ref_seq(:, im), kin.B, 'unit', 'rad');
    joint_q = keni_sol_inverse(T_nom, kin.B, kin.l0, kin.P_m, p_seq);
    q_act = [joint_q(4,1); joint_q(3,2); joint_q(3,3); joint_q(3,4); joint_q(3,5)];

    % 笛卡尔直接逆解对照
    for i = 1:5
        q_direct = norm(T_nom(1:3,4) + T_nom(1:3,1:3)*kin.P_m(:,i) - kin.B(:,i)) - kin.l0(i);
        assert(abs(q_direct - q_act(i)) < 1e-9, '支链%d逆解不一致', i);
    end

    [T_poe, ~] = keni_sol_forward(joint_q, p_seq, 1e-8);
    X0 = trans2dq(T_nom);   % 初值取名义位姿（实际使用方式）
    [T_dq, ~, info] = dq_fks_spr4ups(q_act, kin, X0, 1e-10, 50);
    d_pos = norm(T_poe(1:3,4) - T_dq(1:3,4));
    R_err = T_poe(1:3,1:3).' * T_dq(1:3,1:3);
    d_rot = acos(max(-1, min(1, (trace(R_err)-1)/2)));
    fprintf('pose %d: FKS iter=%d res=%.2e | Δpos=%.3e m, Δrot=%.3e rad\n', ...
        im, info.iter, info.res, d_pos, d_rot);
    assert(info.converged, 'DQ-FKS 未收敛');
    assert(d_pos < 1e-6 && d_rot < 1e-6, 'DQ-FKS 与 POE 正解不一致');
end

%% 2) EMM 解析 vs 有限差分
fprintf('===== 2) EMM 解析 vs 有限差分 =====\n');
im = 1;
kin.d = 0.3e-3;   % 取非零 d 以覆盖非理想约束路径（文献 §2.2）
T_nom = pos2trans(Pos_ref_seq(:, im), kin.B, 'unit', 'rad');
joint_q = keni_sol_inverse(T_nom, kin.B, kin.l0, kin.P_m, p_seq);
q_act = [joint_q(4,1); joint_q(3,2); joint_q(3,3); joint_q(3,4); joint_q(3,5)];
[T0, X0] = dq_fks_spr4ups(q_act, kin, [], 1e-10, 50);

tool.c  = [100; 100; 50] * unit_para;
tool.R0 = eye(3);
tool.o0 = [0; 0; 37.8] * unit_para;

J1 = build_emm_j1_spr4ups(T0, kin);
J3 = build_emm_j3(T0, J1, tool);

h = 1e-6;   % 扰动步长 (m)
err_j1 = zeros(35, 1);
err_j3 = zeros(38, 1);
for k = 1 : 38
    kin_p = kin;
    tool_p = tool;
    if k <= 35
        i = ceil(k/7);
        off = k - 7*(i-1);
        if off <= 3
            kin_p.B(off, i) = kin_p.B(off, i) + h;
        elseif off == 4
            kin_p.l0(i) = kin_p.l0(i) + h;
        elseif i == 1
            % 支链 1 参数块 [Δb1; ΔL1; Δd; Δa1y; Δa1z]（当前坐标系中 Δa1x 由 Δd 替代）
            if off == 5
                kin_p.d = kin_p.d + h;
            elseif off == 6
                kin_p.P_m(2, 1) = kin_p.P_m(2, 1) + h;
            else
                kin_p.P_m(3, 1) = kin_p.P_m(3, 1) + h;
            end
        else
            kin_p.P_m(off-4, i) = kin_p.P_m(off-4, i) + h;
        end
        [T_p, ~] = dq_fks_spr4ups(q_act, kin_p, X0, 1e-10, 50);
        d_o = T_p(1:3,4) - T0(1:3,4);
        R_err = T_p(1:3,1:3) * T0(1:3,1:3).';
        d_th = 0.5 * [R_err(3,2)-R_err(2,3); R_err(1,3)-R_err(3,1); R_err(2,1)-R_err(1,2)];  % 空间系姿态误差一阶形式（ΔR≈[Δθ×]R）
        err_j1(k) = norm([d_o; d_th]/h - J1(:, k));

        m0 = local_feature_points(T0, tool);
        mp = local_feature_points(T_p, tool);
        err_j3(k) = norm((mp - m0)/h - reshape(J3(:, k), 3, 3), 'fro');
    else
        j = k - 35;
        tool_p.c(j) = tool_p.c(j) + h;
        m0 = local_feature_points(T0, tool);
        mp = local_feature_points(T0, tool_p);
        err_j3(k) = norm((mp - m0)/h - reshape(J3(:, k), 3, 3), 'fro');
    end
end
fprintf('J1 各列相对误差 max = %.3e\n', max(err_j1)/max(sqrt(sum(J1.^2,1))));
fprintf('J3 各列相对误差 max = %.3e\n', max(err_j3)/max(sqrt(sum(J3.^2,1))));
assert(max(err_j1) < 1e-3 * max(sqrt(sum(J1.^2,1))), 'J1 与有限差分不符');
assert(max(err_j3) < 1e-3 * max(sqrt(sum(J3.^2,1))), 'J3 与有限差分不符');
%% 3) 非理想约束逆解（文献 §2.2）回环校验
fprintf('===== 3) iks_nonideal 回环校验 =====\n');
for d_test = [0, 0.3e-3, -0.5e-3]
    kin.d = d_test;
    for im = [1, round(size(Pos_ref_seq,2)/2)]
        [q_ik, T_ik] = iks_nonideal(Pos_ref_seq(:, im), kin);
        % d = 0 时应与 pos2trans 完全一致
        if d_test == 0
            T_ref = pos2trans(Pos_ref_seq(:, im), kin.B, 'unit', 'rad');
            assert(norm(T_ik - T_ref) < 1e-12, 'd=0 时 iks_nonideal 与 pos2trans 不一致');
        end
        % 正解回环：由 q 反解位姿应与 IKS 位姿一致
        [T_fk, ~, info] = dq_fks_spr4ups(q_ik, kin, trans2dq(T_ik), 1e-10, 50);
        assert(info.converged, '回环校验：FKS 未收敛');
        fprintf('d=%.1f mm, pose %d: |T_fk - T_ik|=%.3e\n', ...
            d_test*1e3, im, norm(T_fk - T_ik));
        assert(norm(T_fk - T_ik) < 1e-8, 'IKS/FKS 回环不一致（d=%.1f mm）', d_test*1e3);
    end
end
fprintf('全部验证通过。\n');

function m = local_feature_points(T, tool)
    c1 = [0;0;0]; c2 = [tool.c(1);0;0]; c3 = [tool.c(2);tool.c(3);0];
    m = [T(1:3,4) + T(1:3,1:3)*(tool.o0 + tool.R0*c1), ...
         T(1:3,4) + T(1:3,1:3)*(tool.o0 + tool.R0*c2), ...
         T(1:3,4) + T(1:3,1:3)*(tool.o0 + tool.R0*c3)];
end
