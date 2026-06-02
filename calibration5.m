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
fprintf('>>>= calibration5 (RTLS) start (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));


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


err_max = 3e-6;
keni_forward_err_max = 1e-8;
loop_max = 100;
para_err_std = 0.0010;  % 制造误差标准差（m/rad），模拟真实参数与名义参数的偏差

% 外部阻尼参数（类似 LM 的 λ，叠加在 TLS 噪声估计之上）
% lambda_damp = 0      → 纯 TLS（σ_n(A) > σ_{n+1} 时）
% lambda_damp = σ_{n+1}² → 标准 LS
% lambda_damp > σ_{n+1}² → 阻尼 LS（保守步长）
lambda_damp = 0.01;

% ----- input data ------
x_seq = [-100 0 100] * unit_para;
y_seq = [-100 0 100] * unit_para;
z_seq = [-900 -1000] * unit_para;
phi_seq = deg2rad(-120 : 120 : 120); 
theta_seq = deg2rad([0 5 10]);
Pos_ref_seq = zeros(5, 3*3*3*3*2);
for ix = 1:length(x_seq)
    for iy = 1:length(y_seq)
        for iz = 1:length(z_seq)
            for iphi = 1:length(phi_seq)
                for itheta = 1:length(theta_seq)
                    Pos_ref_seq(:, (ix-1)*54+(iy-1)*18+(iz-1)*9+(iphi-1)*3+(itheta)) = [x_seq(ix); y_seq(iy); z_seq(iz); phi_seq(iphi); theta_seq(itheta)];
                end
            end
        end
    end
end

% x_seq = [-100 0 100] * unit_para;
% y_seq = [-100 0 100] * unit_para;
% z_seq = [-800 -850 -900] * unit_para;
% theta_seq = [0 5 10];
% Pos_ref_seq = zeros(5, 3*3*3*3);
% for ix = 1:3
%     for iy = 1:3
%         for iz = 1:3
%             for itheta = 1:3
%                 Pos_ref_seq(:, (ix-1)*27+(iy-1)*9+(iz-1)*3+(itheta)) = [x_seq(ix); y_seq(iy); z_seq(iz); 0; theta_seq(itheta)];
%             end
%         end
%     end
% end
seq_len = length(Pos_ref_seq(1, :));
% ----- end input data ------


%% 标定步骤
[p_seq_nom, xi_seq_nom] = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

% 创建"真实"模型：在名义参数上叠加制造误差
rng(0313);
p_seq_true = p_seq_nom + para_err_std * randn(6, 34);

T_cal_seq = zeros(4, 4, seq_len);
T_measure_seq = zeros(4, 4, seq_len);
joint_seq_iter = zeros(6, 5, seq_len);
err_seq_iter = zeros(6*seq_len, 1);
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

    % 标定初始关节量（名义逆解值）
    joint_seq_iter(:,:,im) = joint_q_ref;

    % 初始位姿误差
    err_seq_iter(6*(im-1)+1 : 6*im) = log_se3(T_measure_seq(:,:,im) / T_cal_seq(:,:,im));
end


%% RTLS 迭代标定
calib_loop = 0;
Jp_bar = zeros(6*seq_len, 112);
err_list = zeros(loop_max, 1);
method_list = cell(loop_max, 1);
lambda_list = zeros(loop_max, 1);

err_cur = norm(err_seq_iter);
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
