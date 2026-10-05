% 参考文献：c13 Kinematic calibration of a 5-axis parallel machining robot
%           based on dimensionless error mapping matrix (Xuan Luo et al., 2021)
% 标定方案：J3 无量纲误差映射矩阵（式(52)~(55)）——3 个靶球特征点的位置
%           误差等效表达末端全位姿误差，并引入辅助工具误差 c，共 38 个
%           误差参数（35 个结构参数 + 3 个工具参数，全部长度量纲）。
%           结构参数含支链 1 非理想约束偏置 d：文献坐标系中 d 沿 e2、替代
%           Δa1y；当前坐标系 x=文献y、y=−文献x、z=文献z，故当前 d 沿 e1、
%           替代 Δa1x，支链 1 参数块为 [b1, L1, d, a1y, a1z]
% 正解方案：对偶四元数 Newton 迭代（文献 §2.2~§2.3 非理想约束模型），
%           见 lib_calib/dq_fks_spr4ups.m
% 辨识流程：文献 §4.1 迭代辨识方案 Step 2~6
% 数据装配：参考 calibration4_exp.m（理论位姿 calib_p/pos_seq.csv，
%           实测靶球点/位姿 calib_p/t1~t3.txt，按序号一一对应）
%
% 说明：文献 Step 1（基于抗扰指标优选测量位姿）属于测量前规划，本脚本
%       处理的是已采集数据，不含该步骤。
clear
path_add();
fprintf('>>>= start (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));

%% 参数集（名义几何，单位 m）
basic_paras = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'm');
unit_para = basic_paras.unit_para;   % mm -> m

kin.B   = basic_paras.B;       % 3×5 基座铰点坐标
kin.P_m = basic_paras.P_m;     % 3×5 动平台铰点坐标
kin.l0  = basic_paras.l0_seq;  % 1×5 零位杆长
kin.d   = 0;                   % 支链 1 非理想约束偏置（文献 §2.2，名义为 0）

% ===================================================================
% 辅助测量工具参数（{at} 工具坐标系 / eeT_at 固定变换）
% -------------------------------------------------------------------
% {at} 建立约定（文献 §3.1/§5.1 的载体）：靶球与文献特征点的对应关系为
%   c_1 = t2（{at} 原点），c_2 = t3，c_3 = t1
% 因此 t2/t3/t1 在 {at} 中的坐标为
%   c_1 = [0;0;0],  c_2 = [c(1);0;0],  c_3 = [c(2);c(3);0]
% eeT_at = [R0, o0]：{at} -> {ee}（动平台表面坐标系）的固定齐次变换。
% 默认按文献式(60)~(63) 由参考位姿的实测数据自动初始化（工具尺寸与
% 安装关系未单独测定时的做法）；若已单独精确测定（如文献 Step 2~3 在
% 零位拟合主轴/参考平面），设 tool_auto_init = false 并手动填写。
% 注意：以下手动填写的值采用 x 轴由 t3 指向 t2 的 {at} 取向（即文献
% 约定绕 z 轴转 180°，故 c(1)<0、R0=I）；feature_points_calc 与
% build_emm_j3 对 c 是线性的，两种取向等价，但同一组 tool 值内部
% 必须自洽。
% ===================================================================
tool_auto_init = false;
tool_ref_pose  = 1;      % 参考位姿序号（建议取接近零位、测量质量好的位姿）
if ~tool_auto_init
    % 原点为 t2（= c_1），x 轴 t3->t2，c_2 = t3 位于 x 轴负向，
    % c_3 = t1 位于 xy 平面第三象限
    tool.c  = [-303.1089; -151.5544; -262.5000] * unit_para;  % 工具名义尺寸 [c1; c2; c3]
    tool.R0 = eye(3);                       % {at} 相对 {ee} 的姿态
    tool.o0 = [151.5544; 87.5000; -37.8] * unit_para;     % {at} 原点在 {ee} 中的坐标
end

% ===================================================================
% Sim Para Config
% ===================================================================
tol_fks     = 1e-10;   % 对偶四元数正解收敛阈值
max_iter_fks = 50;     % 正解最大迭代次数
tol_ident   = 1e-7;    % 参数误差范数收敛阈值（m）
loop_max    = 20;      % 辨识最大迭代次数
lambda_reg  = 0;       % 正则化系数，0 = 普通最小二乘（文献 Step 5）；
                       % >0 时为 Tikhonov 正则化以应对辨识矩阵病态

% ----- input data ------
% 理论位姿：calib_p/pos_seq.csv，每行 [x y z phi theta]，单位 mm / deg
pos_csv = readmatrix(fullfile('calib_p', 'pos_seq.csv'));
Pos_ref_seq = pos_csv.';                                   % 5×n
Pos_ref_seq(1:3, :) = Pos_ref_seq(1:3, :) * unit_para;     % mm -> m
Pos_ref_seq(4:5, :) = deg2rad(Pos_ref_seq(4:5, :));        % deg -> rad

% 测量数据：动平台位姿序列 + 世界系下的原始靶球点（J3 残差用后者）
% T_ref_seq 为 t1~t3 首末行零位复测点构建的位姿（末尾零位对比用）
[~, T_measure_seq, pts_meas_seq, T_ref_seq] = calib_pts2pose_seq('calib_p');
T_measure_seq(1:3, 4, :) = T_measure_seq(1:3, 4, :) * unit_para;  % mm -> m
pts_meas_seq = pts_meas_seq * unit_para;                          % mm -> m
T_ref_seq(1:3, 4, :) = T_ref_seq(1:3, 4, :) * unit_para;          % mm -> m
% 靶球与文献特征点对应关系：c_1 = t2, c_2 = t3, c_3 = t1。
% pts_meas_seq 原始列序为 [t1 t2 t3]，重排为文献序号 [c_1 c_2 c_3]，
% 此后残差堆叠、feature_points_calc、build_emm_j3 均按 c 序号一致处理
pts_meas_seq = pts_meas_seq(:, [2 3 1], :);

seq_len = size(Pos_ref_seq, 2);
assert(size(T_measure_seq, 3) == seq_len, ...
    '测量位姿数(%d)与理论位姿数(%d)不一致', size(T_measure_seq, 3), seq_len);
assert(9*seq_len >= 38, ...
    '测量位姿数(%d)不足：J3 每姿态提供 9 个方程，至少需 5 个位姿', seq_len);
% ----- end input data ------

% 工具参数自动初始化（文献式(60)~(63)：参考位姿下 eeT_at = T_ee⁻¹·T_at）
if tool_auto_init
    tool = init_tool_params(T_measure_seq(:, :, tool_ref_pose), ...
        pts_meas_seq(:, :, tool_ref_pose));
    fprintf('工具参数自动初始化（参考位姿 #%d）：c = [%.4f %.4f %.4f] mm\n', ...
        tool_ref_pose, tool.c(1)/unit_para, tool.c(2)/unit_para, tool.c(3)/unit_para);
end

%% 名义逆解 → 主动关节指令（文献 Step 2：名义驱动输入）
% 非理想约束逆解（文献 §2.2 式(3)~(8)，e2→e1 适配本机构 x 轴约束）：
% 给定 (x,y,z,phi,theta)，由约束 (b1-o)·x = d 确定 R，再计算
%   L1 = |o + R(a1 + d·e1) - b1|,  L_i = |o + R·a_i - b_i| (i=2..5)
% 注：d = 0 时退化为理想约束，与 pos2trans 结果一致
q_act_seq = zeros(5, seq_len);   % 支链长度序列
X_dq_seq  = zeros(8, seq_len);   % 对偶四元数初值/热启动缓存
for im = 1 : seq_len
    [q_act_seq(:, im), T_nom] = iks_nonideal(Pos_ref_seq(:, im), kin);
    X_dq_seq(:, im) = trans2dq(T_nom);   % 首轮正解初值取名义位姿
end

%% 标定前残差（名义参数）
kin_nom = kin;
res_pre = eval_residuals(kin_nom, tool, q_act_seq, X_dq_seq, ...
    T_measure_seq, pts_meas_seq, unit_para, tol_fks, max_iter_fks);
print_residuals('标定前', res_pre);

%% 迭代辨识（文献 §4.1 Step 3~6）
%   Step 3: 当前几何参数 + 名义驱动输入 → FKS 计算理论位姿/特征点坐标
%   Step 4: 实测与理论特征点坐标相减得残差，计算各姿态 EMM 并堆叠
%   Step 5: 最小二乘辨识几何误差，判断范数是否小于阈值
%   Step 6: 未收敛则将误差叠加到当前参数，回到 Step 3
W3 = zeros(9*seq_len, 38);   % 堆叠辨识矩阵（文献式(59)）
X3 = zeros(9*seq_len, 1);    % 堆叠残差向量
err_list = zeros(loop_max+1, 1);

% 误差参数名称（38 维，诊断不可辨识方向用）：
% 支链 1 [Δb1x Δb1y Δb1z ΔL1 Δd Δa1y Δa1z]（当前坐标系中 d 沿 e1 方向，
% 替代 Δa1x），支链 2~5 [Δb_ix Δb_iy Δb_iz ΔL_i Δa_ix Δa_iy Δa_iz]，
% 末尾 [Δc1 Δc2 Δc3]
param_names = cell(38, 1);
for i = 1 : 5
    param_names(7*(i-1) + (1:7)) = cellstr([ ...
        "b_"+i+"x"; "b_"+i+"y"; "b_"+i+"z"; "L_"+i; ...
        "a_"+i+"x"; "a_"+i+"y"; "a_"+i+"z"]);
end
param_names(5:7) = {'d', 'a_1y', 'a_1z'};
param_names(36:38) = {'c_1', 'c_2', 'c_3'};

calib_loop = 0;
rms_prev = res_pre.feat_rmse;
kin_prev = kin;
tool_prev = tool;
while true
    % ---- Step 3 & 4 ----
    for im = 1 : seq_len
        [T_fks, X_dq_seq(:, im), info_fks] = dq_fks_spr4ups(q_act_seq(:, im), kin, ...
            X_dq_seq(:, im), tol_fks, max_iter_fks);
        if ~info_fks.converged
            warning('位姿 #%d 正解未收敛（‖F‖=%.2e），该姿态 EMM 可能不可靠。', ...
                im, info_fks.res);
        end

        % 特征点残差：Δm_i = m_i^实测 - m_i^计算（9×1）
        m_calc = feature_points_calc(T_fks, tool);
        dm = pts_meas_seq(:, :, im) - m_calc;
        X3(9*(im-1)+1 : 9*im) = dm(:);

        % 各姿态无量纲 EMM：J3 = J_at · [J1, 0; 0, I3]
        J1 = build_emm_j1_spr4ups(T_fks, kin);
        W3(9*(im-1)+1 : 9*im, :) = build_emm_j3(T_fks, J1, tool);
    end

    % 发散/收敛判断：残差明显增大（>50%）视为发散；轻微增大视为已到噪声
    % 水平，均回退到上一轮参数并终止
    rms_cur = rms(X3);
    if calib_loop > 0 && rms_cur > rms_prev  % rms_prev只有calib_loop=0时是距离平均值，后续都是上一轮的rms_cur
        kin = kin_prev;
        tool = tool_prev;
        calib_loop = calib_loop - 1;
        if rms_cur > rms_prev * 1.5
            warning('辨识发散：残差 %.4f -> %.4f mm，回退到上一轮参数并终止。', ...
                rms_prev/unit_para, rms_cur/unit_para);
        else
            fprintf('残差已到噪声水平（%.4f mm），终止迭代。\n', rms_prev/unit_para);
        end
        break;
    end
    err_list(calib_loop+1) = rms_cur / unit_para;

    % ---- Step 5：最小二乘辨识（截断 SVD，滤除不可辨识方向）----
    % 无量纲 EMM 缓解了混合量纲导致的病态，但不可辨识方向（如测量位姿
    % 激励不足）仍需通过奇异值截断或 Tikhonov 正则化处理
    if lambda_reg > 0
        eps_ident = (W3.'*W3 + lambda_reg*eye(38)) \ (W3.'*X3);
        r_eff = NaN;
        sv = [NaN; NaN];
    else
        [U_s, S_s, V_s] = svd(W3, 'econ');
        sv = diag(S_s);
        sv_tol = 1e-8 * sv(1);
        idx = sv > sv_tol;
        r_eff = sum(idx);
        eps_ident = V_s(:, idx) * ((U_s(:, idx).' * X3) ./ sv(idx));
        if r_eff < 38
            % 最弱可辨识方向的主要参数构成（可用于指导增补测量位姿）
            [~, order] = sort(abs(V_s(:, end)), 'descend');
            fprintf('  最弱可辨识方向主要参数: %s\n', ...
                strjoin(param_names(order(1:5)), ', '));
        end
    end

    calib_loop = calib_loop + 1;
    fprintf(['loop = %d, 特征点残差 RMSE = %.4f mm, ‖Δε‖ = %.3e m, ' ...
        'rank(W3) = %d/38, σmax/σmin = %.2e\n'], calib_loop, rms_cur/unit_para, ...
        norm(eps_ident), r_eff, sv(1)/sv(end));

    if norm(eps_ident) < tol_ident || calib_loop >= loop_max
        break;
    end

    % ---- Step 6：误差叠加到当前几何参数 ----
    kin_prev = kin;
    tool_prev = tool;
    rms_prev = rms_cur;
    % 支链 1 参数块 [Δb1(3); ΔL1; Δd; Δa1y; Δa1z]（当前坐标系中 Δa1x 由 Δd 替代）
    kin.B(:, 1)   = kin.B(:, 1)   + eps_ident(1:3);
    kin.l0(1)     = kin.l0(1)     + eps_ident(4);
    kin.d         = kin.d         + eps_ident(5);
    kin.P_m(2, 1) = kin.P_m(2, 1) + eps_ident(6);
    kin.P_m(3, 1) = kin.P_m(3, 1) + eps_ident(7);
    for i = 2 : 5
        cidx = 7*(i-1) + (1:7);             % 支链 i 参数块 [Δb; ΔL; Δa]
        kin.B(:, i)   = kin.B(:, i)   + eps_ident(cidx(1:3));
        kin.l0(i)     = kin.l0(i)     + eps_ident(cidx(4));
        kin.P_m(:, i) = kin.P_m(:, i) + eps_ident(cidx(5:7));
    end
    tool.c = tool.c + eps_ident(36:38);
end

%% 标定后残差
res_post = eval_residuals(kin, tool, q_act_seq, X_dq_seq, ...
    T_measure_seq, pts_meas_seq, unit_para, tol_fks, max_iter_fks);
print_residuals('标定后', res_post);

%% 标定前后残差对比图
fig_val = figure('Color', [1 1 1]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(1:seq_len, res_pre.dpos, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0); hold on;
plot(1:seq_len, res_post.dpos, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
yline(mean(res_pre.dpos), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
yline(mean(res_post.dpos), '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 0.8);
ylabel('Δd (mm)', 'FontSize', 12, 'FontName', 'Times New Roman');
legend({'标定前', '标定后'}, 'FontSize', 11, 'FontName', '微软雅黑', 'Location', 'best');
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'LineWidth', 1.0);
grid on; box on;

nexttile;
plot(1:seq_len, res_pre.dang, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0); hold on;
plot(1:seq_len, res_post.dang, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
yline(mean(res_pre.dang), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
yline(mean(res_post.dang), '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 0.8);
ylabel('Δθ (°)', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'LineWidth', 1.0);
grid on; box on;

nexttile;
plot(1:seq_len, res_pre.dfeat, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0); hold on;
plot(1:seq_len, res_post.dfeat, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
yline(mean(res_pre.dfeat), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
yline(mean(res_post.dfeat), '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 0.8);
xlabel('Pos labels', 'FontSize', 12, 'FontName', '微软雅黑');
ylabel('Δp_{feat} (mm)', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'LineWidth', 1.0);
grid on; box on;

fig = figure('Color', [1 1 1]);
plot(0:calib_loop, err_list(1:calib_loop+1), 'linewidth', 1.5)
set(gca, 'YScale', 'log');
grid on
set(gca, 'linewidth', 1.5, 'fontsize', 15, 'fontname', 'Times New Roman');
set(gcf, 'unit', 'centimeters', 'position', [10 10 14 8]);
xlabel('迭代次数', 'FontSize', 14, 'FontName', '微软雅黑', 'FontWeight', 'bold');
ylabel('特征点残差 RMSE (mm)', 'FontSize', 14, 'FontName', '微软雅黑', 'FontWeight', 'bold');

%% 零位位姿（基于标定后参数）
% 零位 = 各主动关节输入为零（支链长 = l0），由标定后参数经对偶四元数
% 正解求得；初值取参考位姿（接近零位）正解的热启动结果
[T_zero, ~, info_zero] = dq_fks_spr4ups(zeros(5, 1), kin, ...
    X_dq_seq(:, tool_ref_pose), tol_fks, max_iter_fks);
if ~info_zero.converged
    warning('零位正解未收敛（‖F‖=%.2e），导出的零位位姿可能不可靠。', info_zero.res);
end
% 与 pos2trans/calib_pts2pose_seq 一致的姿态约定：由 R 的 z 轴列恢复 phi/theta
pos_zero = [T_zero(1:3, 4) / unit_para; ...
            atan2d(T_zero(2, 3), T_zero(1, 3)); ...
            acosd(max(-1, min(1, T_zero(3, 3))))];   % [x y z] mm, [phi theta] deg
fprintf(['零位位姿（标定后）: x=%.4f  y=%.4f  z=%.4f mm,  ' ...
    'phi=%.4f°  theta=%.4f°\n'], pos_zero);

%% 标定参数导出为 CSV（单位: mm）
% 注意：本脚本输出的是 c13 式笛卡尔参数集（B, l0, P_m, d, c），与
% calibrated_params_exp.csv 的 POE 导出格式不同，如需接入现有控制器
% 需经 parameterize() 转换（P_m 的 xy 分量需与 r1/r2、limb_dir 协调；
% 且控制器的正逆解需采用相同的非理想约束模型，d 才有效）
calib_csv = 'calibrated_params_c13_exp.csv';
fid = fopen(calib_csv, 'w');
if fid < 0
    warning('无法写入 %s（文件可能被其他程序占用），本次跳过 CSV 导出。', calib_csv);
else
fprintf(fid, '# SPR-4UPS calibrated kinematic parameters (c13 dimensionless EMM, J3)\n');
fprintf(fid, '# Units: mm (length), deg (angle)\n');
fprintf(fid, '# Columns: param_name, value[, value...]\n');
fprintf(fid, '# zero_pose: platform pose at q=0 from calibrated parameters, [x,y,z] mm, [phi,theta] deg\n');
for i = 1:5
    fprintf(fid, 'B_%d,%.12f,%.12f,%.12f\n', i, ...
        kin.B(1,i)/unit_para, kin.B(2,i)/unit_para, kin.B(3,i)/unit_para);
end
for i = 1:5
    fprintf(fid, 'l0_%d,%.12f\n', i, kin.l0(i)/unit_para);
end
for i = 1:5
    fprintf(fid, 'Pm_%d,%.12f,%.12f,%.12f\n', i, ...
        kin.P_m(1,i)/unit_para, kin.P_m(2,i)/unit_para, kin.P_m(3,i)/unit_para);
end
fprintf(fid, 'd,%.12f\n', kin.d/unit_para);
fprintf(fid, 'tool_c,%.12f,%.12f,%.12f\n', ...
    tool.c(1)/unit_para, tool.c(2)/unit_para, tool.c(3)/unit_para);
fprintf(fid, 'zero_pose,%.12f,%.12f,%.12f,%.12f,%.12f\n', ...
    pos_zero(1), pos_zero(2), pos_zero(3), pos_zero(4), pos_zero(5));
fclose(fid);
fprintf('标定参数已导出至 %s\n', calib_csv);
end

%% 零位实测位姿（t1~t3 首末行复测点）与标定后零位正解对比
% 首末行为零位参考点的两次复测，T_ref_seq 由三点构建（世界系，与标定中
% T_measure_seq 同一约定，即视测量世界系与基座系一致）；T_zero 为标定后
% 参数 q=0 的正解零位位姿
fprintf('零位对比（首末行实测复测点 vs 标定后正解零位）——\n');
fprintf('  正解零位: [%.4f %.4f %.4f] mm, phi=%.4f° theta=%.4f°\n', pos_zero);
for k = 1 : 2
    T_ref = T_ref_seq(:, :, k);
    dt = (T_ref(1:3, 4) - T_zero(1:3, 4)) / unit_para;    % mm
    dz = acosd(max(-1, min(1, dot(T_ref(1:3, 3), T_zero(1:3, 3)))));  % 主轴向夹角
    R_err = T_ref(1:3, 1:3).' * T_zero(1:3, 1:3);
    dR = rad2deg(acos(max(-1, min(1, (trace(R_err) - 1) / 2))));      % 全姿态误差
    phi_ref = atan2d(T_ref(2, 3), T_ref(1, 3));
    theta_ref = acosd(max(-1, min(1, T_ref(3, 3))));
    if k == 1, tag = '首行'; else, tag = '末行'; end
    fprintf(['  %s: 实测 [%.4f %.4f %.4f] mm, phi=%.4f° theta=%.4f° | ' ...
        'Δpos = %.4f mm ([%+.4f %+.4f %+.4f]), Δ主轴向 = %.4f°, Δ全姿态 = %.4f°\n'], ...
        tag, T_ref(1:3, 4)/unit_para, phi_ref, theta_ref, ...
        norm(dt), dt, dz, dR);
end
% 首末行互差（零位复测重复性）
dt_rep = (T_ref_seq(1:3, 4, 2) - T_ref_seq(1:3, 4, 1)) / unit_para;
dz_rep = acosd(max(-1, min(1, dot(T_ref_seq(1:3, 3, 1), T_ref_seq(1:3, 3, 2)))));
fprintf('  首末行互差（重复性）: Δpos = %.4f mm, Δ主轴向 = %.4f°\n', ...
    norm(dt_rep), dz_rep);

fprintf('>>>= done (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));


% ======================================================================
% 局部函数
% ======================================================================
function m_calc = feature_points_calc(T, tool)
% 由当前位姿计算 3 个靶球特征点在基座系下的坐标（文献式(51)）
%   m_i = o + R·(o0 + R0·c_i),  c_1=[0;0;0], c_2=[c1;0;0], c_3=[c2;c3;0]
    o = T(1:3, 4);
    R = T(1:3, 1:3);
    c1 = [0; 0; 0];
    c2 = [tool.c(1); 0; 0];
    c3 = [tool.c(2); tool.c(3); 0];
    m_calc = [o + R*(tool.o0 + tool.R0*c1), ...
              o + R*(tool.o0 + tool.R0*c2), ...
              o + R*(tool.o0 + tool.R0*c3)];
end


function res = eval_residuals(kin, tool, q_act_seq, X_dq_seq, ...
    T_measure_seq, pts_meas_seq, unit_para, tol_fks, max_iter_fks)
% 逐姿态评估残差：TCP 位置误差、主轴姿态误差、靶球特征点误差
    seq_len = size(q_act_seq, 2);
    res.dpos  = zeros(seq_len, 1);   % TCP 位置误差 (mm)
    res.dang  = zeros(seq_len, 1);   % 姿态误差 (°)
    res.dfeat = zeros(seq_len, 1);   % 特征点位置误差均值 (mm)
    for im = 1 : seq_len
        [T_cal, ~, info_fks] = dq_fks_spr4ups(q_act_seq(:, im), kin, ...
            X_dq_seq(:, im), tol_fks, max_iter_fks);
        if ~info_fks.converged
            warning('残差评估：位姿 #%d 正解未收敛（‖F‖=%.2e）。', im, info_fks.res);
        end

        res.dpos(im) = norm(T_measure_seq(1:3,4,im) - T_cal(1:3,4)) / unit_para;
        z_meas = T_measure_seq(1:3, 3, im);
        z_cal  = T_cal(1:3, 3);
        res.dang(im) = rad2deg(acos(max(min(dot(z_meas, z_cal), 1), -1)));

        m_err = pts_meas_seq(:, :, im) - feature_points_calc(T_cal, tool);
        res.dfeat(im) = mean(sqrt(sum(m_err.^2, 1))) / unit_para;  % 靶球空间距离误差平均值，三维坐标误差sum(x, 1)，对x的列向量求和
    end
    res.feat_rmse = rms(res.dfeat) * unit_para;   % 内部单位 m
end


function print_residuals(tag, res)
    fprintf(['%s ——\n pos:  mean=%.4f  max=%.4f  rmse=%.4f mm, \n' ...
        ' ang:  mean=%.4f  max=%.4f  rmse=%.4f°,\n' ...
        ' feat: mean=%.4f  max=%.4f  rmse=%.4f mm\n'], tag, ...
        mean(res.dpos), max(res.dpos), rms(res.dpos), ...
        mean(res.dang), max(res.dang), rms(res.dang), ...
        mean(res.dfeat), max(res.dfeat), rms(res.dfeat));
end


function tool = init_tool_params(T_surf, qpts)
% 由参考位姿的实测数据初始化辅助工具参数（文献式(60)~(63)）
%   T_surf — 参考位姿下动平台表面坐标系的实测位姿（4×4）
%   qpts   — 参考位姿下 3 个特征点在世界系下的实测坐标（3×3），列序已按
%            文献序号排列 [c_1 c_2 c_3] = [t2 t3 t1]（在主脚本中完成重排）
%   {at}：原点 c_1，x 轴 c_1→c_2，c_3 在 xy 平面内（文献约定，c(1)>0）；
%   eeT_at = T_surf⁻¹·T_at
%   注：此取向与主脚本中手动填写的 tool 值（x 轴 t3→t2，c(1)<0）相差
%   绕 z 轴 180° 旋转，二者等价但数值不同，混用前需统一到同一取向。
    q1 = qpts(:, 1);
    q2 = qpts(:, 2);
    q3 = qpts(:, 3);
    x_at = (q2 - q1) / norm(q2 - q1);
    z_at = cross(x_at, q3 - q1);
    z_at = z_at / norm(z_at);
    y_at = cross(z_at, x_at);
    T_at = [[x_at, y_at, z_at], q1; 0 0 0 1];

    eeT_at = T_surf \ T_at;
    tool.R0 = eeT_at(1:3, 1:3);
    tool.o0 = eeT_at(1:3, 4);
    tool.c  = [norm(q2 - q1);
               dot(q3 - q1, x_at);
               dot(q3 - q1, y_at)];
end
