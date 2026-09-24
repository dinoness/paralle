function [J1, J_ee, J_1mat] = build_emm_j1_spr4ups(T, kin)
%BUILD_EMM_J1_SPR4UPS 传统误差映射矩阵 J1（c13 式(44)/(45) 在 SPR-4UPS 上的对应）
%   位姿误差 δ1 = [Δo; Δθ]（Δθ 为等效转轴误差矢量，ΔR ≈ [Δθ×]R）与结构
%   误差参数 ε 之间的映射：
%       δ1 = J1·ε,   J1 = J_ee \ J_1mat
%   ε 排序（每条支链 7 个，共 35 个，全部为长度量纲；与文献一致，支链 1
%   的非理想约束偏置 Δd 替代一个动平台铰点分量——文献坐标系中 d 沿 e2、
%   替代 Δa1y；当前坐标系 x=文献y、y=−文献x、z=文献z，故当前 d 沿 e1、
%   替代 Δa1x，保留的铰点误差为 Δa1y 和 Δa1z，文献 §2.2/§3.1）：
%       ε_1 = [Δb_1(3); ΔL_1(1); Δd; Δa1y; Δa1z]
%       ε_i = [Δb_i(3); ΔL_i(1); Δa_i(3)],  i = 2..5
%   J_ee / J_1mat 的行构成（J_1mat 取运动学方程对参数偏导的负号约定，
%   与本文件历史版本一致）：
%     第 1 行 — 支链 1 长方程 b1 + L1 s1 = o + R(a1 + d·e1)（文献式(43)，
%       连接点 a1* = a1 + d·e1）微分后投影到支链方向 s1：
%         J_ee(1,:)  = [s1ᵀ,  (R a1* × s1)ᵀ]
%         J_1mat(1, 支链1块) = [s1ᵀ,  1,  -s1ᵀRe1,  -s1ᵀRe2,  -s1ᵀRe3]
%           （第 5 列对应 Δd，系数 -s1ᵀRe1；第 6 列对应 Δa1y，系数
%             -s1ᵀRe2；第 7 列对应 Δa1z，系数 -s1ᵀRe3）
%     第 i 行（i=2..5）— 支链长方程 b_i + L_i s_i = o + R a_i：
%         J_ee(i,:)  = [s_iᵀ,  (R a_i × s_i)ᵀ]
%         J_1mat(i, 支链i块) = [s_iᵀ,  1,  -s_iᵀ R]
%     第 6 行 — 非理想 SPR 约束 (b_1 - o)·(R e1) = d（文献式(43)，
%       当前坐标系中 d 沿 e1 方向；d=0 退化为理想约束）微分：
%         J_ee(6,:)  = [-(R e1)ᵀ,  ((R e1)×(b1-o))ᵀ]
%         J_1mat(6, 支链1块) = [-(R e1)ᵀ,  0,  1,  0,  0]
%           （第 5 列对应 Δd，系数 +1）
%
%   输入：
%     T   — 4×4 当前动平台位姿（FKS 计算值）
%     kin — 几何参数结构体（.B, .l0, .P_m, .d，见 dq_fks_spr4ups；
%           .d 缺省按 0 处理）
%   输出：
%     J1     — 6×35 传统 EMM（混合量纲）
%     J_ee   — 6×6 位姿误差系数阵
%     J_1mat — 6×35 结构误差系数阵

    o = T(1:3, 4);
    R = T(1:3, 1:3);
    d = 0;
    if isfield(kin, 'd')
        d = kin.d;
    end
    e1 = [1; 0; 0];
    e2 = [0; 1; 0];
    e3 = [0; 0; 1];

    J_ee = zeros(6, 6);
    J_1mat = zeros(6, 35);

    % 支链行（i = 1..5）
    for i = 1 : 5
        a = kin.P_m(:, i);
        if i == 1
            a = a + d*e1;               % 支链 1 实际连接点 a1*（文献式(43)）
        end
        w = o + R*a - kin.B(:, i);
        s = w / norm(w);                    % 支链单位方向矢量
        Ra_cross_s = cross(R*a, s);

        J_ee(i, :) = [s.', Ra_cross_s.'];

        cidx = 7*(i-1) + (1:7);             % 支链 i 在 ε 中的列块
        if i == 1
            % ε_1 = [Δb1(3); ΔL1; Δd; Δa1y; Δa1z]
            J_1mat(i, cidx) = [s.', 1, -s.'*R*e1, -s.'*R*e2, -s.'*R*e3];
        else
            J_1mat(i, cidx) = [s.', 1, -s.'*R];
        end
    end

    % 非理想 SPR 约束行：(b1 - o)·(R e1) = d
    u = R * e1;
    J_ee(6, :) = [-u.', cross(u, kin.B(:,1) - o).'];
    J_1mat(6, 1:7) = [-u.', 0, 1, 0, 0];

    if rcond(J_ee) < 1e-12
        warning('build_emm_j1_spr4ups:IllConditioned', ...
            'J_ee 接近奇异（机构可能接近奇异位形），J1 可能不可靠。');
    end
    J1 = J_ee \ J_1mat;
end
