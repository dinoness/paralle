% 参考文献：c02 并联机构的运动学误差建模及参数可辨识性分析_孔令雨
% L-M求解
% 实验版本：理论位姿来自 calib_p/pos_seq.csv，测量位姿由
% lib_calib/calib_pts2pose_seq.m 从 calib_p/t1~t3.txt 实测点重构，
% 两者按序号一一对应
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
% -----end-struct-parameter------

% ===================================================================
% Sim Para Config
% ===================================================================
err_max = 2e-5;
keni_forward_err_max = 1e-8;
loop_max = 100;

% ----- input data ------
% 理论位姿：calib_p/pos_seq.csv，每行 [x y z phi theta]，单位 mm / deg
pos_csv = readmatrix(fullfile('calib_p', 'pos_seq.csv'));
Pos_ref_seq = pos_csv.';                                   % 5×n
Pos_ref_seq(1:3, :) = Pos_ref_seq(1:3, :) * unit_para;     % mm -> m
Pos_ref_seq(4:5, :) = deg2rad(Pos_ref_seq(4:5, :));        % deg -> rad

% 测量位姿：由实测靶球点构建的坐标系直接生成齐次变换矩阵序列，
% 4×4×n，平移单位 mm，与理论位姿按序号一一对应
[~, T_measure_seq] = calib_pts2pose_seq('calib_p');
T_measure_seq(1:3, 4, :) = T_measure_seq(1:3, 4, :) * unit_para;   % mm -> m

seq_len = size(Pos_ref_seq, 2);
assert(size(T_measure_seq, 3) == seq_len, ...
    '测量位姿数(%d)与理论位姿数(%d)不一致', size(T_measure_seq, 3), seq_len);
% ----- end input data ------

%% 标定步骤
[p_seq_nom, xi_seq_nom] = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

T_cal_seq = zeros(4, 4, seq_len);
joint_seq_iter = zeros(6, 5, seq_len);
err_seq_iter = zeros(6*seq_len, 1);

p_seq_iter = p_seq_nom;  % 结构参数序列
xi_seq_iter = xi_seq_nom;

%% 实验数据装配
%   理论位姿 pos_seq.csv → 名义逆解 → 关节量（指令值）
%   测量位姿由 calib_pts2pose_seq 以变换矩阵形式直接给出
for im = 1 : seq_len
    T_cal_seq(:,:,im) = pos2trans(Pos_ref_seq(:, im), B, 'unit', 'rad');

    % 名义模型逆解得到所有关节量（含主动+被动）
    joint_seq_iter(:,:,im) = keni_sol_inverse(T_cal_seq(:,:,im), B, l0_seq, P_m, p_seq_nom);

    % 初始位姿误差
    err_seq_iter(6*(im-1)+1 : 6*im) = log_se3(T_measure_seq(:,:,im) / T_cal_seq(:,:,im));
end

%% 引入真实位姿，判断误差是否小于容差，是则结束，否则进入迭代
calib_loop = 0;
Jp_bar = zeros(6*seq_len, 112);
err_list = zeros(loop_max+1, 1);
lambda = 1e-1;

err_cur = norm(err_seq_iter);
err_list(1) = err_cur;

%% 标定前残差（名义参数 p_seq_nom）
dpos_pre = zeros(seq_len, 1);   % 位置误差 (mm)
dang_pre = zeros(seq_len, 1);   % 姿态误差 (°)
for im = 1 : seq_len
    [T_cal_pre, ~] = keni_sol_forward(joint_seq_iter(:,:,im), p_seq_nom, keni_forward_err_max);
    pos_meas = trans2pos(T_measure_seq(:,:,im));
    pos_cal  = trans2pos(T_cal_pre);
    dpos_pre(im) = norm(pos_meas(1:3) - pos_cal(1:3)) / unit_para;
    z_meas = T_measure_seq(1:3, 3, im);
    z_cal  = T_cal_pre(1:3, 3);
    dang_pre(im) = rad2deg(acos(max(min(dot(z_meas, z_cal), 1), -1)));
end
fprintf("标定前 ——\n pos:  mean=%.4f  max=%.4f  rmse=%.4f mm, \n ang:  mean=%.4f  max=%.4f  rmse=%.4f°\n", ...
    mean(dpos_pre), max(dpos_pre), rms(dpos_pre), mean(dang_pre), max(dang_pre), rms(dang_pre));

while err_cur > err_max
    [U, N, V_prep, M] = calib_iter_row_space_matrix(xi_seq_iter);
    [Lambda, Ap] = calib_iter_restore_matrix(p_seq_iter);
    % 更新运动学参数
    for im = 1 : seq_len
        Jp_bar(6*(im-1)+1 : 6*(im-1) + 6, :) = calib_iter_matrix2(joint_seq_iter(:,:,im), p_seq_iter, xi_seq_iter, U, V_prep);
    end


    pk = (Jp_bar' * Jp_bar + lambda * eye(size(Jp_bar, 2))) \ (Jp_bar' * err_seq_iter);

    delta_p_seq = restore_full_param_increment(pk, U, V_prep, Lambda, Ap);  % 依据式(3-28)恢复为原始运动学参数
    p_seq_vec = p_seq_iter(:) + delta_p_seq(:);  % A(:)矩阵按列排列成列向量
    p_seq_iter = reshape(p_seq_vec, size(p_seq_iter, 1), size(p_seq_iter, 2));  % 将向量重排为矩阵

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

%% 标定后残差（标定后参数 p_seq_iter）
dpos_post = zeros(seq_len, 1);   % 位置误差 (mm)
dang_post = zeros(seq_len, 1);   % 姿态误差 (°)
for im = 1 : seq_len
    [T_cal_post, ~] = keni_sol_forward(joint_seq_iter(:,:,im), p_seq_iter, keni_forward_err_max);
    pos_meas = trans2pos(T_measure_seq(:,:,im));
    pos_cal  = trans2pos(T_cal_post);
    dpos_post(im) = norm(pos_meas(1:3) - pos_cal(1:3)) / unit_para;
    z_meas = T_measure_seq(1:3, 3, im);
    z_cal  = T_cal_post(1:3, 3);
    dang_post(im) = rad2deg(acos(max(min(dot(z_meas, z_cal), 1), -1)));
end
fprintf("标定后 ——\n pos:  mean=%.4f  max=%.4f  rmse=%.4f mm, \n ang:  mean=%.4f  max=%.4f  rmse=%.4f°\n", ...
    mean(dpos_post), max(dpos_post), rms(dpos_post), mean(dang_post), max(dang_post), rms(dang_post));

%% 标定前后残差对比图
fig_val = figure('Color', [1 1 1]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(1:seq_len, dpos_pre, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0); hold on;
plot(1:seq_len, dpos_post, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
yline(mean(dpos_pre), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
yline(mean(dpos_post), '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 0.8);
ylabel('Δd (mm)', 'FontSize', 12, 'FontName', 'Times New Roman');
legend({'标定前', '标定后'}, 'FontSize', 11, 'FontName', '微软雅黑', 'Location', 'best');
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'LineWidth', 1.0);
grid on; box on;

nexttile;
plot(1:seq_len, dang_pre, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0); hold on;
plot(1:seq_len, dang_post, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0);
yline(mean(dang_pre), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
yline(mean(dang_post), '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 0.8);
xlabel('Pos labels', 'FontSize', 12, 'FontName', '微软雅黑');
ylabel('Δθ (°)', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'LineWidth', 1.0);
grid on; box on;

fig = figure('Color', [1 1 1]);
plot(0:calib_loop, err_list(1:calib_loop+1), 'linewidth', 1.5)
set(gca, 'YScale', 'log');  % 对数坐标轴
grid on
set(gca, 'linewidth', 1.5, 'fontsize', 15, 'fontname', 'Times New Roman');
set(gcf, 'unit', 'centimeters', 'position', [10 10 14 8]);
xlabel('迭代次数', 'FontSize', 14, 'FontName', '微软雅黑', 'FontWeight', 'bold', 'Color', 'black');
ylabel('残余误差', 'FontSize', 14, 'FontName', '微软雅黑', 'FontWeight', 'bold', 'Color', 'black');

[B_out, r1_out, r2_out, l0_out, P_m_out, limb_out, joint_u_angle_tilt_out] = ...
    deparameterize(p_seq_iter, P_m, limb_dir);

%% 标定参数导出为 CSV（单位: mm, °）
calib_csv = 'calibrated_params_exp.csv';
fid = fopen(calib_csv, 'w');
fprintf(fid, '# SPR-4UPS calibrated kinematic parameters (experimental data)\n');
fprintf(fid, '# Units: mm (length), deg (angle)\n');
fprintf(fid, '# Columns: param_name, value[, value...]\n');
fprintf(fid, '# Rows: one parameter per line, label then values\n');
fprintf(fid, 'r1,%.12f\n', r1_out / unit_para);
fprintf(fid, 'r2,%.12f\n', r2_out / unit_para);
fprintf(fid, 'joint_u_angle_tilt,%.12f\n', rad2deg(joint_u_angle_tilt_out));
for i = 1:5
    fprintf(fid, 'l0_%d,%.12f\n', i, l0_out(i) / unit_para);
end
for i = 1:5
    fprintf(fid, 'B_%d,%.12f,%.12f,%.12f\n', i, B_out(1,i) / unit_para, B_out(2,i) / unit_para, B_out(3,i) / unit_para);
end
for i = 1:5
    fprintf(fid, 'Pm_%d,%.12f,%.12f,%.12f\n', i, P_m_out(1,i) / unit_para, P_m_out(2,i) / unit_para, P_m_out(3,i) / unit_para);
end
for i = 1:5
    fprintf(fid, 'limb_dir_%d,%.12f,%.12f\n', i, rad2deg(limb_out(i,1)), rad2deg(limb_out(i,2)));
end
fclose(fid);
fprintf('标定参数已导出至 %s\n', calib_csv);

fprintf('>>>= done (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));
