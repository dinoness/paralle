function [q_act, T] = iks_nonideal(pos, kin)
%IKS_NONIDEAL 非理想约束逆解（文献 c13 §2.2 式(3)~(8)，e2→e1 适配本机构）
%   给定指令位姿 pos = [x; y; z; phi; theta]（m / rad，姿态约定与
%   pos2trans 一致），由非理想约束 (b1 - o)·x = d 确定姿态矩阵 R：
%     z 轴由 (phi, theta) 确定；x ⊥ z 且 x·(b1-o) = d（d=0 时退化为
%     理想约束，结果与 pos2trans 一致）；y = z×x。
%   支链长（文献式(8)，支链 1 连接点 a1* = a1 + d·e1）：
%     L1 = |o + R(a1 + d·e1) - b1|,  L_i = |o + R·a_i - b_i|,  i = 2..5
%
%   输入：
%     pos   — 5×1 指令位姿 [x; y; z; phi; theta]，长度 m，角度 rad
%     kin   — 几何参数结构体（.B, .l0, .P_m, .d，见 dq_fks_spr4ups）
%   输出：
%     q_act — 5×1 主动关节量（q_i = L_i - l0_i）
%     T     — 4×4 动平台位姿（基座系）
    t = pos(1:3);
    phi = pos(4);
    theta = pos(5);
    z = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

    w = kin.B(:, 1) - t;               % b1 - o
    s0 = cross(w, z);
    rho = norm(s0);                    % w 在 ⊥z 平面内的投影长度
    s0 = s0 / rho;                     % 理想约束下的 x 轴（同 pos2trans）
    t0 = cross(z, s0);                 % ⊥z 平面内另一单位方向（t0·w = rho）
    sd = kin.d / rho;                  % sinψ：x 相对理想方向的偏角
    assert(abs(sd) < 1, '非理想约束逆解：|d|=%.3e m 超出可行范围（rho=%.3e m）', ...
        kin.d, rho);
    x = sqrt(1 - sd^2)*s0 + sd*t0;     % 取与理想约束连续的分支（cosψ>0）
    y = cross(z, x);
    R = [x, y, z];
    T = [R, t; 0 0 0 1];

    q_act = zeros(5, 1);
    for i = 1 : 5
        a = kin.P_m(:, i);
        if i == 1
            a = a + kin.d*[1; 0; 0];
        end
        q_act(i) = norm(t + R*a - kin.B(:, i)) - kin.l0(i);
    end
end
