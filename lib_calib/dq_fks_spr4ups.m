function [T, X, info] = dq_fks_spr4ups(q_act, kin, X0, tol, max_iter)
%DQ_FKS_SPR4UPS 基于单位对偶四元数的 SPR-4UPS 运动学正解（Newton 迭代）
%   参考文献 c13 §2.2~§2.3（式(22)~(32)）：以单位对偶四元数
%       X = [qv(3); q0; pv(3); p0]  (8×1)
%   为广义坐标，正解转化为求解关于 X 的二次方程组 F(X)=0：
%       f_1 = |o + R·(a1 + d·e1) - b1|^2 - L1^2 = 0     （支链 1 连接点带
%                                                         非理想偏置 d，式(22)）
%       f_i = |o + R·a_i - b_i|^2 - L_i^2 = 0,  i = 2..5
%       f_6 = (b_1 - o)·(R·e1) = d                       （非理想 SPR 约束，
%                                                         式(24) 取 e2→e1 适配
%                                                         本机构 x 轴约束；
%                                                         d=0 退化为理想约束，
%                                                         与 pos2trans 一致）
%       f_7 = q·q - 1 = 0,  f_8 = q·p = 0                 （单位对偶四元数约束）
%   文献的快速迭代式(32) X←½X+D⁻¹C 与一般 Newton 步 X←X-D⁻¹F 在 F 为纯
%   二次型时代数等价；本实现直接由对偶四元数运算给出解析的 F 与 Jacobian
%   D=∂F/∂X，采用带回退的阻尼 Newton 步，二次收敛且对初值更稳健。
%
%   输入：
%     q_act    — 5×1 主动关节量（各支链移动副位移）
%     kin      — 几何参数结构体：
%                  .B    3×5 基座铰点坐标（基座系）
%                  .l0   1×5 各支链零位杆长
%                  .P_m  3×5 动平台铰点坐标（平台系）
%                  .d    (可选) 支链 1 非理想约束偏置（m），缺省 0
%     X0       — (可选) 8×1 迭代初值（对偶四元数），缺省取平台零位
%     tol      — (可选) 收敛阈值，缺省 1e-10
%     max_iter — (可选) 最大迭代次数，缺省 50
%   输出：
%     T    — 4×4 动平台位姿（基座系）
%     X    — 8×1 收敛的对偶四元数（可作下一组输入的热启动初值）
%     info — 结构体：.iter 迭代次数, .res 收敛时 ‖F‖, .converged 是否收敛
%
%   位姿与对偶四元数的关系（文献式(20)(21)，平移未取 1/2 因子）：
%     o = Im(p·q̄) = q0·pv - p0·qv + qv×pv
%     R·a = q·a·q̄ = (q0²-qv·qv)·a + 2(qv·a)qv + 2q0(qv×a)

    if nargin < 3 || isempty(X0)
        X0 = [0; 0; 0; 1; 0; 0; 0; 0];  % 单位四元数 + 零平移
    end
    if nargin < 4 || isempty(tol)
        tol = 1e-10;
    end
    if nargin < 5 || isempty(max_iter)
        max_iter = 50;
    end

    L = kin.l0(:) + q_act(:);   % 各支链实际杆长
    d = 0;                      % 非理想约束偏置（文献 §2.2）
    if isfield(kin, 'd')
        d = kin.d;
    end

    X = X0;
    converged = false;
    [F, D] = fk_equations(X, kin.B, kin.P_m, L, d);
    res = norm(F);
    for it = 1 : max_iter
        if res < tol
            converged = true;
            break;
        end
        if rcond(D) < 1e-14
            warning('dq_fks_spr4ups:SingularD', ...
                'D 接近奇异（机构可能接近奇异位形），迭代终止。');
            break;
        end
        % 阻尼 Newton：全步长残差增大则回退半步，提高初值较差时的稳健性
        dX = -D \ F;
        alpha = 1;
        for ls = 1 : 10
            X_new = X + alpha * dX;
            [F_new, D_new] = fk_equations(X_new, kin.B, kin.P_m, L, d);
            res_new = norm(F_new);
            if res_new < res
                break;
            end
            alpha = alpha / 2;
        end
        X = X_new;  F = F_new;  D = D_new;  res = res_new;
    end

    % 由收敛的对偶四元数重建位姿（q 归一化）
    qv = X(1:3);  q0 = X(4);  pv = X(5:7);  p0 = X(8);
    nq = hypot(norm(qv), q0);
    qv = qv / nq;  q0 = q0 / nq;
    o = q0*pv - p0*qv + cross(qv, pv);
    R = rot_apply(qv, q0, eye(3));
    T = [R, o; 0 0 0 1];

    info.iter = it;
    info.res = res;
    info.converged = converged;
end


function [F, D] = fk_equations(X, B, P_m, L, d)
% 正解方程组 F(X) 及其解析 Jacobian D = ∂F/∂X
    qv = X(1:3);  q0 = X(4);  pv = X(5:7);  p0 = X(8);

    o = q0*pv - p0*qv + cross(qv, pv);
    % o 对各分量的偏导
    do_dqv = -p0*eye(3) - skew(pv);
    do_dq0 = pv;
    do_dpv = q0*eye(3) + skew(qv);
    do_dp0 = -qv;

    F = zeros(8, 1);
    D = zeros(8, 8);

    % f_1..f_5：支链长方程（支链 1 连接点 a1* = a1 + d·e1，文献式(22)(23)）
    for i = 1 : 5
        a = P_m(:, i);
        if i == 1
            a = a + d*[1; 0; 0];
        end
        [Ra, dRa_dqv, dRa_dq0] = rot_apply_grad(qv, q0, a);  % Ra为旋转变换后的a
        w = o + Ra - B(:, i);
        F(i) = w.'*w - L(i)^2;  % F1-F5的定义式支链的长度误差
        dw = [do_dqv + dRa_dqv, do_dq0 + dRa_dq0, do_dpv, do_dp0];  % 3×8
        D(i, :) = 2*w.'*dw;  % 对于F1-F5的表达式对X求导
    end

    % f_6：非理想 SPR 约束 (b1 - o)·(R·e1) = d（文献式(24)，e2→e1）
    e1 = [1; 0; 0];
    [u, du_dqv, du_dq0] = rot_apply_grad(qv, q0, e1);
    g = B(:, 1) - o;
    F(6) = g.'*u - d;
    D(6, :) = [-u.'*do_dqv + g.'*du_dqv, ...
               -u.'*do_dq0 + g.'*du_dq0, ...
               -u.'*do_dpv, -u.'*do_dp0];

    % f_7, f_8：单位对偶四元数约束
    F(7) = qv.'*qv + q0^2 - 1;
    D(7, :) = [2*qv.', 2*q0, 0, 0, 0, 0];
    F(8) = qv.'*pv + q0*p0;
    D(8, :) = [pv.', p0, qv.', q0];
end


function Ra = rot_apply(qv, q0, A)
% 对 A 的各列施加旋转 R·a = q·a·q̄（q 未归一化时结果带 |q|² 比例）
    Ra = (q0^2 - qv.'*qv).*A + 2*qv*(qv.'*A) + 2*q0*skew(qv)*A;
end


function [Ra, dRa_dqv, dRa_dq0] = rot_apply_grad(qv, q0, a)
% R·a 及其对 qv、q0 的解析偏导
    Ra = rot_apply(qv, q0, a);
    dRa_dqv = 2*(qv.'*a)*eye(3) + 2*qv*a.' - 2*a*qv.' - 2*q0*skew(a);
    dRa_dq0 = 2*q0*a + 2*cross(qv, a);
end
