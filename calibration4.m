% 参考文献：c02 并联机构的运动学误差建模及参数可辨识性分析_孔令雨
% L-M求解
clear
path_add();
fprintf('>>>= start (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));


%% 参数集
% -----struct-parameter------
basic_paras = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'm');  % 单位意思是程序中用到的单位
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
% -----end-struct-parameter------

err_max = 2e-5;
keni_forward_err_max = 1e-8;
loop_max = 50;
para_err_std = 0.0005;  % 制造误差标准差（m/rad），模拟真实参数与名义参数的偏差
noise_measure_std = 2e-5;  % 跟踪仪误差5um/m


% ----- input data ------
[Pos_ref_seq, Pos_valid_seq] = calib_seq_generate(unit_para);
seq_len = size(Pos_ref_seq, 2);
valid_len = size(Pos_valid_seq, 2);
% ----- end input data ------

%% 标定步骤
[p_seq_nom, xi_seq_nom] = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

% 创建"真实"模型：在名义参数上叠加制造误差
rng(0313);
% p_seq_true = p_seq_nom + para_err_std * randn(6, 34);
p_seq_true = p_seq_nom + para_err_std * rand(6, 34);

T_cal_seq = zeros(4, 4, seq_len);
T_measure_seq = zeros(4, 4, seq_len);
joint_seq_iter = zeros(6, 5, seq_len);
err_seq_iter = zeros(6*seq_len, 1);

T_measure_valid_seq = zeros(4, 4, valid_len);
joint_valid_seq = zeros(6, 5, valid_len);
p_seq_iter = p_seq_nom;  % 结构参数序列
xi_seq_iter = xi_seq_nom;

%% 仿真测量数据生成
%   步骤1：名义逆解 → 主动关节量（"指令值"）
%   步骤2：真实模型正解 → 实际输出位姿（"测量值"）
for im = 1 : seq_len
    T_cal_seq(:,:,im) = pos2trans(Pos_ref_seq(:, im), B,'unit','rad');

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

%% 引入真实位姿，判断误差是否小于容差，是则结束，否则进入迭代
calib_loop = 0;
Jp_bar = zeros(6*seq_len, 112);
err_list = zeros(loop_max+1,1);
lambda = 1e-1;

err_cur = norm(err_seq_iter);
err_list(1) = err_cur;

%% 标定前验证集残差（名义参数 p_seq_nom）
err_valid_pre = zeros(6*valid_len, 1);
dpos_pre = zeros(valid_len, 1);   % 位置误差 (mm)
dang_pre = zeros(valid_len, 1);   % 姿态误差 (°)
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
    % size(U{1})
    % size(V_prep{1})
    % size(V{1})
    [Lambda, Ap] = calib_iter_restore_matrix(p_seq_iter);
    % 更新运动学参数
    for im = 1 : seq_len
        % Jp_bar(6*(im-1)+1 : 6*(im-1) + 6, :) = calib_iter_matrix(joint_seq_iter(:,:,im), p_seq_iter);
        Jp_bar(6*(im-1)+1 : 6*(im-1) + 6, :) = calib_iter_matrix2(joint_seq_iter(:,:,im), p_seq_iter, xi_seq_iter, U, V_prep);
    end
    

    pk = (Jp_bar' * Jp_bar + lambda * eye(size(Jp_bar, 2))) \ (Jp_bar' * err_seq_iter);
    
    delta_p_seq = restore_full_param_increment(pk, U, V_prep, Lambda, Ap);  % 依据式(3-28)恢复为原始运动学参数
    p_seq_vec = p_seq_iter(:) + delta_p_seq(:);  % A(:)矩阵按列排列成列向量
    p_seq_iter = reshape(p_seq_vec, size(p_seq_iter, 1), size(p_seq_iter, 2));  % 将向量重排为矩阵
    % p_seq_vec = p_seq_iter(:) + 0.001 * pinv(Jp_bar' * Jp_bar) * Jp_bar' * err_seq_iter;
    

    xi_seq_iter = rebuild_xi_seq_from_p(p_seq_iter);  % 更新零位全局旋量坐标，供下一轮迭代使用
    
    for im = 1 : seq_len
        % 新运动学参数下的正解
        [T_cal_seq(:,:,im), joint_seq_iter(:,:,im)] = keni_sol_forward(joint_seq_iter(:,:,im), p_seq_iter, keni_forward_err_max);
        % 计算位姿误差
        err_seq_iter(6*(im-1)+1 : 6*im) = log_se3(T_measure_seq(:,:,im) / T_cal_seq(:,:,im));
    end

    
    % 阻尼系数迭代
    err_new = norm(err_seq_iter);
    delta_err = err_new - err_cur;  % 实际下降
    delta_qk = -1 *( (Jp_bar'*err_seq_iter)'*pk + 0.5*pk'*(Jp_bar'*Jp_bar)*pk );  % 期望下降
    eta = delta_err / delta_qk;
            
    % 系数可调整
    if eta > 0.75
        lambda = lambda * 0.5;
    elseif eta < 0.25
        lambda = lambda * 3;
    end
    

    calib_loop = calib_loop + 1;
    err_cur = err_new;
    err_list(calib_loop+1) = err_cur;
    
    if(rem(calib_loop, 10) == 0)
        fprintf("loop = %d\n", calib_loop);
    end
    if calib_loop >= loop_max
        break;
    end
    
end

%% 标定后验证集残差（标定后参数 p_seq_iter）
err_valid_post = zeros(6*valid_len, 1);
dpos_post = zeros(valid_len, 1);   % 位置误差 (mm)
dang_post = zeros(valid_len, 1);   % 姿态误差 (°)
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
ylabel('Δd (mm)', 'FontSize', 12, 'FontName', 'Times New Roman');
legend({'标定前', '标定后'}, 'FontSize', 11, 'FontName', '微软雅黑', 'Location', 'best');
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'LineWidth', 1.0);
grid on; box on;

nexttile;
plot(1:valid_len, dang_pre, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0); hold on;
plot(1:valid_len, dang_post, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
yline(mean(dang_pre), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
yline(mean(dang_post), '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 0.8);
xlabel('Valid pos labels', 'FontSize', 12, 'FontName', '微软雅黑');
ylabel('Δθ (°)', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'LineWidth', 1.0);
grid on; box on;

fig = figure('Color', [1 1 1]);
plot(0:calib_loop, err_list(1:calib_loop+1),'linewidth',1.5)
set(gca, 'YScale', 'log');  % 对数坐标轴
grid on
set(gca,'linewidth',1.5,'fontsize',15,'fontname','Times New Roman');
set(gcf,'unit','centimeters','position',[10 10 14 8]);
xlabel('迭代次数', 'FontSize', 14, 'FontName', '微软雅黑', 'FontWeight', 'bold', 'Color', 'black');
ylabel('残余误差', 'FontSize', 14, 'FontName', '微软雅黑', 'FontWeight', 'bold', 'Color', 'black');

fprintf('>>>= done (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));