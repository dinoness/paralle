% 参考文献：c02 并联机构的运动学误差建模及参数可辨识性分析_孔令雨 第四章
% 使用正则化总体最小二乘法 (Regularized Total Least Squares) 进行标定求解
% 与 calibration4.m（LM 阻尼最小二乘）对应
%
% 求解策略（组合正则化）：
%   x = V * diag(s_i / (s_i² - σ_{n+1}² + λ_ext)) * U' * b
% - λ_ext = 0 且 TLS 适定时：纯 TLS（式 4-18），对 Jacobian 误差进行补偿
% - λ_ext > σ_{n+1}²：阻尼最小二乘（类似 LM），保证步长稳健
% - λ_ext 通过增益比 η 自适应调整（与 calibration4.m 策略一致）
clear
path_add();
fprintf('>>>= start (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));


%% 参数集
basic_paras = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'm');
l_max = basic_paras.l_max;
l_min = basic_paras.l_min;
R1 = basic_paras.R1;
R2 = basic_paras.R2;
H = basic_paras.H;
r1 = basic_paras.r1;
r2 = basic_paras.r2;
h = basic_paras.h;
L_tool = basic_paras.L_tool;
limb_dir = basic_paras.limb_dir;
B = basic_paras.B;
P_m = basic_paras.P_m;
l0_seq = basic_paras.l0_seq;
joint_u_angle_tilt = basic_paras.joint_u_angle_tilt;
unit_para = basic_paras.unit_para;

% 标定板参数（标定靶球中心对于动平台的坐标）
Rc = 150 * unit_para;
hc = -50 * unit_para;
calib_dir = deg2rad([-90 30 150]);
Pc = [Rc*cos(calib_dir);
      Rc*sin(calib_dir);
      hc hc hc];

% 标定靶球三点定义的测量坐标系M → 动平台坐标系的固定变换
% T_platform_M: X_platform = T_platform_M * X_M
T_platform2calib = three_pts2trans(Pc(:,1), Pc(:,2), Pc(:,3));

err_max = 2e-5;
keni_forward_err_max = 1e-8;
loop_max = 100;
para_err_std = 0.0005;  % 制造误差标准差（m/rad），模拟真实参数与名义参数的偏差
noise_measure_std = 5e-6;  % 跟踪仪误差5um/m

% 外部阻尼参数（类似 LM 的 λ，叠加在 TLS 噪声估计之上）
% lambda_damp = 0      → 纯 TLS（σ_n(A) > σ_{n+1} 时）
% lambda_damp = σ_{n+1}² → 标准 LS
% lambda_damp > σ_{n+1}² → 阻尼 LS（保守步长）
lambda_damp = 0.01;

% ----- input data ------
[Pos_ref_seq, Pos_valid_seq] = calib_seq_generate(unit_para);
seq_len = size(Pos_ref_seq, 2);
valid_len = size(Pos_valid_seq, 2);


%% 标定步骤
[p_seq_nom, xi_seq_nom] = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

% 创建"真实"模型：在名义参数上叠加制造误差
rng(0313);
p_seq_true = p_seq_nom + para_err_std * rand(6, 34);

T_cal_seq = zeros(4, 4, seq_len);
T_measure_seq = zeros(4, 4, seq_len);
joint_seq_iter = zeros(6, 5, seq_len);
err_seq_iter = zeros(6*seq_len, 1);

T_measure_valid_seq = zeros(4, 4, valid_len);
joint_valid_seq = zeros(6, 5, valid_len);
p_seq_iter = p_seq_nom;
xi_seq_iter = xi_seq_nom;

%% 仿真测量数据生成
%   步骤1：名义逆解 → 主动关节量（"指令值"）
%   步骤2：真实模型正解 → 实际输出位姿（"测量值"）
for im = 1 : seq_len
    T_cal_seq(:,:,im) = pos2trans(Pos_ref_seq(:, im), B);

    % 步骤1：名义模型逆解得到所有关节量（含主动+被动）
    joint_q_ref = keni_sol_inverse(T_cal_seq(:,:,im), B, l0_seq, P_m, p_seq_nom);

    % 步骤2：固定主动关节（移动副），用真实模型正解得到实际位姿
    [T_measure_seq(:,:,im), ~] = keni_sol_forward(joint_q_ref, p_seq_true, keni_forward_err_max);

    % 标定靶球在世界系下的坐标 + 高斯噪声 → 重构含噪位姿
    Pc_world = T_measure_seq(:,:,im) * [Pc; 1, 1, 1];          % 4×3，末行全1
    Pc_world_noisy = Pc_world(1:3, :) + noise_measure_std * randn(3, 3);
    T_M_world_noisy = three_pts2trans(Pc_world_noisy(:,1), Pc_world_noisy(:,2), Pc_world_noisy(:,3));
    T_measure_seq(:,:,im) = T_M_world_noisy / T_platform2calib;

    % 标定初始关节量（名义逆解值）
    joint_seq_iter(:,:,im) = joint_q_ref;

    % 初始位姿误差
    err_seq_iter(6*(im-1)+1 : 6*im) = log_se3(T_measure_seq(:,:,im) / T_cal_seq(:,:,im));
end

%% 验证集生成：与标定集相同的噪声添加流程
for iv = 1 : valid_len
    T_cal_valid = pos2trans(Pos_valid_seq(:, iv), B, 'unit', 'rad');
    joint_q_valid = keni_sol_inverse(T_cal_valid, B, l0_seq, P_m, p_seq_nom);
    [T_measure_valid_seq(:,:,iv), joint_valid_seq(:,:,iv)] = keni_sol_forward(joint_q_valid, p_seq_true, keni_forward_err_max);

    Pc_world = T_measure_valid_seq(:,:,iv) * [Pc; 1, 1, 1];
    Pc_world_noisy = Pc_world(1:3, :) + noise_measure_std * randn(3, 3);
    T_M_world_noisy = three_pts2trans(Pc_world_noisy(:,1), Pc_world_noisy(:,2), Pc_world_noisy(:,3));
    T_measure_valid_seq(:,:,iv) = T_M_world_noisy / T_platform2calib;
end


%% RTLS 迭代标定
calib_loop = 0;
Jp_bar = zeros(6*seq_len, 112);
err_list = zeros(loop_max, 1);
method_list = cell(loop_max, 1);
lambda_list = zeros(loop_max, 1);

err_cur = norm(err_seq_iter);

%% 标定前验证集残差（名义参数 p_seq_nom）
err_valid_pre = zeros(6*valid_len, 1);
dpos_pre = zeros(valid_len, 1);
dang_pre = zeros(valid_len, 1);
for iv = 1 : valid_len
    [T_cal_valid, ~] = keni_sol_forward(joint_valid_seq(:,:,iv), p_seq_nom, keni_forward_err_max);
    err_valid_pre(6*(iv-1)+1 : 6*iv) = log_se3(T_measure_valid_seq(:,:,iv) / T_cal_valid);
    pos_meas = trans2pos(T_measure_valid_seq(:,:,iv));
    pos_cal  = trans2pos(T_cal_valid);
    dpos_pre(iv) = norm(pos_meas(1:3) - pos_cal(1:3)) / unit_para;
    z_meas = T_measure_valid_seq(1:3, 3, iv);
    z_cal  = T_cal_valid(1:3, 3);
    dang_pre(iv) = rad2deg(acos(max(min(dot(z_meas, z_cal), 1), -1)));
end
fprintf("验证集标定前 ——\n pos:  mean=%.4f  max=%.4f  rmse=%.4f mm, \n ang:  mean=%.4f  max=%.4f  rmse=%.4f°\n", ...
    mean(dpos_pre), max(dpos_pre), rms(dpos_pre), mean(dang_pre), max(dang_pre), rms(dang_pre));

while err_cur > err_max
    [U, N, V_prep, M] = calib_iter_row_space_matrix(xi_seq_iter);
    [Lambda, Ap] = calib_iter_restore_matrix(p_seq_iter);

    for im = 1 : seq_len
        Jp_bar(6*(im-1)+1 : 6*(im-1) + 6, :) = calib_iter_matrix2(joint_seq_iter(:,:,im), p_seq_iter, xi_seq_iter, U, V_prep);
    end

    % 组合 TLS + 外部阻尼求解
    [pk, tls_info] = solve_tls(Jp_bar, err_seq_iter, lambda_damp);

    % 保存回退状态（正解失败时回退）
    p_seq_prev = p_seq_iter;
    xi_seq_prev = xi_seq_iter;
    err_prev = err_seq_iter;

    % 直接应用全步长（与 calibration4.m 一致，不进行线搜索）
    delta_p_seq = restore_full_param_increment(pk, U, V_prep, Lambda, Ap);
    p_seq_vec = p_seq_iter(:) + delta_p_seq(:);
    p_seq_iter = reshape(p_seq_vec, size(p_seq_iter, 1), size(p_seq_iter, 2));
    xi_seq_iter = rebuild_xi_seq_from_p(p_seq_iter);

    % 正解验证
    fk_success = true;
    T_cal_try = zeros(4, 4, seq_len);
    joint_q_try = zeros(6, 5, seq_len);
    err_try = zeros(6*seq_len, 1);
    for im = 1 : seq_len
        try
            [T_cal_try(:,:,im), joint_q_try(:,:,im)] = keni_sol_forward(joint_seq_iter(:,:,im), p_seq_iter, keni_forward_err_max);
            err_try(6*(im-1)+1 : 6*im, 1) = log_se3(T_measure_seq(:,:,im) / T_cal_try(:,:,im));
        catch
            fk_success = false;
            break;
        end
    end

    if ~fk_success || any(isnan(err_try)) || any(isinf(err_try))
        % 正解失败 → 增大阻尼，回退参数，重试
        lambda_damp = max(lambda_damp * 10, 1e-3);
        p_seq_iter = p_seq_prev;
        xi_seq_iter = xi_seq_prev;
        err_cur = norm(err_prev);
        calib_loop = calib_loop + 1;
        err_list(calib_loop) = err_cur;
        method_list{calib_loop} = 'FAIL';
        lambda_list(calib_loop) = lambda_damp;
        if calib_loop > loop_max; break; end
        continue;
    end

    % 增益比自适应阻尼（与 calibration4.m LM 策略一致）
    err_new = norm(err_try);
    delta_err = err_new - err_cur;
    delta_qk = -1 * ((Jp_bar' * err_seq_iter)' * pk + 0.5 * pk' * (Jp_bar' * Jp_bar) * pk);

    if abs(delta_qk) > 1e-15
        eta = delta_err / delta_qk;
    else
        eta = 0;
    end

    if eta > 0.75
        lambda_damp = lambda_damp * 0.5;
    elseif eta < 0.25
        lambda_damp = max(lambda_damp * 3, 1e-8);
    end

    % 始终接受步长（与 calibration4.m 一致）
    T_cal_seq = T_cal_try;
    joint_seq_iter = joint_q_try;
    err_seq_iter = err_try;
    err_cur = err_new;

    calib_loop = calib_loop + 1;
    err_list(calib_loop) = err_cur;
    method_list{calib_loop} = tls_info.method;
    lambda_list(calib_loop) = lambda_damp;

    if calib_loop > loop_max
        break;
    end

    if rem(calib_loop, 10) == 0
        fprintf("loop = %d, err = %.2e, lambda = %.1e, %s (rank=%d, gap=%.1e)\n", ...
                calib_loop, err_cur, lambda_damp, tls_info.method, tls_info.eff_rank, tls_info.gap);
    end
end

fprintf('>>>>> RTLS 标定完成: %d 次迭代, 最终误差 = %.3e <<<<<\n', calib_loop, err_cur);

%% 标定后验证集残差（标定后参数 p_seq_iter）
err_valid_post = zeros(6*valid_len, 1);
dpos_post = zeros(valid_len, 1);
dang_post = zeros(valid_len, 1);
for iv = 1 : valid_len
    [T_cal_valid, ~] = keni_sol_forward(joint_valid_seq(:,:,iv), p_seq_iter, keni_forward_err_max);
    err_valid_post(6*(iv-1)+1 : 6*iv) = log_se3(T_measure_valid_seq(:,:,iv) / T_cal_valid);
    pos_meas = trans2pos(T_measure_valid_seq(:,:,iv));
    pos_cal  = trans2pos(T_cal_valid);
    dpos_post(iv) = norm(pos_meas(1:3) - pos_cal(1:3)) / unit_para;
    z_meas = T_measure_valid_seq(1:3, 3, iv);
    z_cal  = T_cal_valid(1:3, 3);
    dang_post(iv) = rad2deg(acos(max(min(dot(z_meas, z_cal), 1), -1)));
end
fprintf("验证集标定后 ——\n pos:  mean=%.4f  max=%.4f  rmse=%.4f mm, \n ang:  mean=%.4f  max=%.4f  rmse=%.4f°\n", ...
    mean(dpos_post), max(dpos_post), rms(dpos_post), mean(dang_post), max(dang_post), rms(dang_post));

%% 标定前后验证集残差对比图
fig_val = figure('Color', [1 1 1]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(1:valid_len, dpos_pre, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0); hold on;
plot(1:valid_len, dpos_post, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
yline(mean(dpos_pre), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
yline(mean(dpos_post), '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 0.8);
ylabel('\Deltad (mm)', 'FontSize', 12, 'FontName', 'Times New Roman');
legend({'标定前', '标定后'}, 'FontSize', 11, 'FontName', '微软雅黑', 'Location', 'best');
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'LineWidth', 1.0);
grid on; box on;

nexttile;
plot(1:valid_len, dang_pre, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0); hold on;
plot(1:valid_len, dang_post, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
yline(mean(dang_pre), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
yline(mean(dang_post), '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 0.8);
xlabel('Valid pos labels', 'FontSize', 12, 'FontName', '微软雅黑');
ylabel('\Delta\theta (°)', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'LineWidth', 1.0);
grid on; box on;


%% 绘图
tls_idx = find(strcmp(method_list(1:calib_loop), 'TLS'));
rtls_idx = find(strcmp(method_list(1:calib_loop), 'RTLS'));
dls_idx = find(strcmp(method_list(1:calib_loop), 'DLS'));

fig = figure('Color', [1 1 1]);

subplot(1, 3, 1);
semilogy(0:calib_loop-1, err_list(1:calib_loop), 'b-', 'linewidth', 1.5);
set(gca, 'linewidth', 1.5, 'fontsize', 12, 'fontname', 'Times New Roman');
xlabel('迭代次数', 'FontSize', 14, 'FontName', '微软雅黑', 'FontWeight', 'bold');
ylabel('残余误差', 'FontSize', 14, 'FontName', '微软雅黑', 'FontWeight', 'bold');
title('标定收敛曲线', 'FontSize', 14, 'FontName', '微软雅黑');
grid on;

subplot(1, 3, 2);
hold on;
if ~isempty(tls_idx)
    stem(tls_idx, ones(size(tls_idx))*1.2, 'b', 'linewidth', 1.5, 'Marker', 'none');
end
if ~isempty(rtls_idx)
    stem(rtls_idx, ones(size(rtls_idx))*1.1, 'g', 'linewidth', 1.5, 'Marker', 'none');
end
if ~isempty(dls_idx)
    stem(dls_idx, ones(size(dls_idx))*1.0, 'r', 'linewidth', 1.5, 'Marker', 'none');
end
hold off;
set(gca, 'linewidth', 1.5, 'fontsize', 12, 'fontname', 'Times New Roman', ...
         'ytick', [1.0 1.1 1.2], 'yticklabel', {'DLS', 'RTLS', 'TLS'});
xlabel('迭代次数', 'FontSize', 14, 'FontName', '微软雅黑', 'FontWeight', 'bold');
title('每轮求解方法', 'FontSize', 14, 'FontName', '微软雅黑');
ylim([0.85 1.35]);
grid on;

subplot(1, 3, 3);
semilogy(1:calib_loop, lambda_list(1:calib_loop), 'r-o', 'linewidth', 1.2, 'MarkerSize', 4);
set(gca, 'linewidth', 1.5, 'fontsize', 12, 'fontname', 'Times New Roman');
xlabel('迭代次数', 'FontSize', 14, 'FontName', '微软雅黑', 'FontWeight', 'bold');
ylabel('阻尼 \lambda', 'FontSize', 14, 'FontName', '微软雅黑', 'FontWeight', 'bold');
title('阻尼自适应', 'FontSize', 14, 'FontName', '微软雅黑');
grid on;

set(gcf, 'unit', 'centimeters', 'position', [10 10 30 8]);

fprintf('>>>= done (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));