clear
fprintf('>>>= Joint_angle_search Start (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));
workspace_plot_flag = 0;  % 1开启工作区间绘制
fig_rotation_show = 0;  % 1开启展示旋转
gif_generate_flag = 0;  % 1为开启录制功能，运行一次程序后记得改文件名


%--------parameter3--------
% 依照 workspace_discrete_v3.m 修改参数配置方式
basic_paras = basic_read('parameters.xlsx', 'column', 'M', 'unit', 'm');

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
unit_para = basic_paras.unit_para;

pos_plant = [0; 0; -800]*unit_para;  % 后面作图用，不参与空间搜索
alpha_plant_deg = 0;
beta_plant_deg = 0;
gamma_plant_deg = 0;
alpha_plant = alpha_plant_deg / 180 * pi;  % 绕 x
beta_plant = beta_plant_deg / 180 * pi;  % 绕 y
gamma_plant = gamma_plant_deg / 180 * pi;  % 绕 z

% static plant (已从 basic_paras 获取 B)
% B 的结构在 basic_read 中已经构建好

% Position and Posture of Move Plant
Rx = [1                0                 0;
      0 cos(alpha_plant) -sin(alpha_plant);
      0 sin(alpha_plant)  cos(alpha_plant)];

Ry = [ cos(beta_plant) 0 sin(beta_plant);
                     0 1               0;
      -sin(beta_plant) 0 cos(beta_plant)];

Rz = [cos(gamma_plant) -sin(gamma_plant) 0;
      sin(gamma_plant)  cos(gamma_plant) 0;
                     0                0  1];

R_plant = Rz * Ry * Rx;

% move plant parameter (已从 basic_paras 获取 P_m)
P_v = zeros(3, 5);  % 只变换了方向，没变换起点
P = zeros(3, 5);    % 末端点坐标
for i = 1 : 5
    P_v(:, i) = R_plant * P_m(:, i);
    P(:, i) = P_v(:, i) + pos_plant;
end

% ball screw dir vector
% 依照 workspace_discrete_v3.m：
% 第一行 phi=35°（球铰链轴线与动平台z轴正向夹角）
% 第二行 theta=limb_dir(:,2)'（方位角与动平台关节一致）
% 之后 ball_vector = -ball_vector，使轴线方向与支链方向同向
% 注：phi=35° + 取反 等价于 原代码中的 phi=145°（不取反）
ball_screw_dir_angle = [deg2rad([35 35 35 35 35]);
                        limb_dir(:, 2)'];
ball_vector = zeros(3, 5);
ball_vector_world = zeros(3, 5);

static_joint_dir_angle_deg = [ 0  0  0  0  0;
                              -90  30 150 -135 -45];
static_joint_dir_angle = static_joint_dir_angle_deg / 180 * pi;
static_joint_vector = zeros(3, 5);

for i_ball = 1 : 5
    ball_vector(1, i_ball) = sin(ball_screw_dir_angle(1, i_ball)) * cos(ball_screw_dir_angle(2, i_ball));
    ball_vector(2, i_ball) = sin(ball_screw_dir_angle(1, i_ball)) * sin(ball_screw_dir_angle(2, i_ball));
    ball_vector(3, i_ball) = cos(ball_screw_dir_angle(1, i_ball));

    static_joint_vector(1, i_ball) = sin(static_joint_dir_angle(1, i_ball)) * cos(static_joint_dir_angle(2, i_ball));
    static_joint_vector(2, i_ball) = sin(static_joint_dir_angle(1, i_ball)) * sin(static_joint_dir_angle(2, i_ball));
    static_joint_vector(3, i_ball) = cos(static_joint_dir_angle(1, i_ball));
end
ball_vector = -1 * ball_vector;  % 与v3一致，让关节方向向量和支链方向同向
ball_vector_world = R_plant * ball_vector;

% -----end-parameter3------



% ===== 阶段1：粗搜索收集各支链方向向量分布 =====
% 使用较粗网格快速收集方向向量，用于后续各支链独立 phi 优化
fprintf('\n>>> Phase 1: Collecting limb direction distributions...\n');
seq_x_coarse = (-400 : 20 : 400)*unit_para;
seq_y_coarse = (-400 : 20 : 400)*unit_para;
seq_z_coarse = (-1200 : 10 : -600)*unit_para;

s_limb_seq = cell(1, 5);
for i = 1 : 5
    s_limb_seq{i} = [];
end

for ix = 1 : length(seq_x_coarse)
    for iy = 1 : length(seq_y_coarse)
        for iz = 1 : length(seq_z_coarse)
            vt = [seq_x_coarse(ix); seq_y_coarse(iy); seq_z_coarse(iz)];
            pos_flag = 0;
            s_limb = zeros(3, 5);
            
            for j = 1 : 5
                vAa = vt + P_v(:, j) - B(:, j);
                len_vAa = norm(vAa);
                s_limb(:, j) = vAa / len_vAa;

                if (len_vAa >= l_min)&&(len_vAa <= l_max)
                    pos_flag = pos_flag + 1;
                end
            end

            if pos_flag == 5
                for i_limb = 1 : 5
                    s_limb_move_in = R_plant' * s_limb(:, i_limb);
                    s_limb_seq{i_limb} = [s_limb_seq{i_limb} s_limb_move_in];
                end
            end
        end
    end
end

for i_limb = 1 : 5
    fprintf('  Limb %d: collected %d direction vectors.\n', i_limb, size(s_limb_seq{i_limb}, 2));
end


% ===== 阶段2：计算各支链最优倾斜角度 phi =====
% 原有方案（v3默认值）：所有支链 phi = 35°
phi_orig = deg2rad([35 35 35 35 35]);

% 新方案：基于方向向量分布，各支链独立计算最优 phi
% 优化目标：使取反后的球铰链轴线方向（-ball_vector）尽可能接近方向向量中心
% 取反后的 ball_vector = -[sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)]
% f(phi) = dot(v_center, -ball_vector) = -sin(phi)*A - cos(phi)*B
% 其中 A = v_x*cos(theta) + v_y*sin(theta), B = v_z
% 最优 phi = atan2(-A, -B)
phi_opt = zeros(1, 5);
for i_limb = 1 : 5
    s_seq = s_limb_seq{i_limb};
    if isempty(s_seq)
        phi_opt(i_limb) = phi_orig(i_limb);
        continue;
    end
    
    % 计算方向向量的平均方向（中心方向）
    v_mean = mean(s_seq, 2);
    v_center = v_mean / norm(v_mean);
    
    % 该支链的 theta（方位角）固定为 limb_dir(i_limb, 2)
    theta = limb_dir(i_limb, 2);
    
    % 最优 phi 解析解
    A = v_center(1) * cos(theta) + v_center(2) * sin(theta);
    B_val = v_center(3);
    phi_calc = atan2(-A, -B_val);
    
    % 将 phi 限制在合理范围 [20°, 60°]（物理可实现范围，且与v3的35°相近）
    if phi_calc < deg2rad(20)
        phi_calc = deg2rad(20);
    elseif phi_calc > deg2rad(60)
        phi_calc = deg2rad(60);
    end
    
    phi_opt(i_limb) = phi_calc;
end

fprintf('\n>>> Phase 2: Optimal phi calculation done.\n');
fprintf('  Original scheme: phi = [%.2f, %.2f, %.2f, %.2f, %.2f] deg\n', rad2deg(phi_orig));
fprintf('  Optimized scheme: phi = [%.2f, %.2f, %.2f, %.2f, %.2f] deg\n', rad2deg(phi_opt));


% ===== 阶段3：评估并对比不同方案（以圆柱体积为优化目标） =====
% 使用中等密度网格评估各方案，并检索最大内接圆柱体积
fprintf('\n>>> Phase 3: Evaluating schemes with cylinder extraction...\n');

% 评估网格（与v3一致）
seq_x_eval = (-400 : 10 : 400)*unit_para;
seq_y_eval = (-400 : 10 : 400)*unit_para;
seq_z_eval = (-1100 : 5 : -700)*unit_para;
angle_limit_deg = 30;  % 球铰链最大允许角度，与v3一致

theta_vec = limb_dir(:, 2)';

% 评估原有方案
[V_orig, R_orig, H_orig, x0_orig, y0_orig, num_orig, max_angles_orig, angle_data_orig] = ...
    eval_scheme(phi_orig, theta_vec, seq_x_eval, seq_y_eval, seq_z_eval, P_v, B, l_min, l_max, R_plant, angle_limit_deg);

% 评估优化方案
[V_opt, R_opt, H_opt, x0_opt, y0_opt, num_opt, max_angles_opt, angle_data_opt] = ...
    eval_scheme(phi_opt, theta_vec, seq_x_eval, seq_y_eval, seq_z_eval, P_v, B, l_min, l_max, R_plant, angle_limit_deg);

% 输出评估结果
fprintf('\n========== Scheme Evaluation (Cylinder Objective) ==========\n');
fprintf('[Original] phi = [35.00, 35.00, 35.00, 35.00, 35.00] deg\n');
fprintf('  Valid workspace points: %d\n', num_orig);
fprintf('  Cylinder: V=%.8f m^3, R=%.6f m, H=%.6f m, center=(%.4f, %.4f) m\n', V_orig, R_orig, H_orig, x0_orig, y0_orig);
fprintf('  Max ball-joint angles:  %.2f, %.2f, %.2f, %.2f, %.2f deg\n', max_angles_orig);

fprintf('\n[Optimized] phi = [%.2f, %.2f, %.2f, %.2f, %.2f] deg\n', rad2deg(phi_opt));
fprintf('  Valid workspace points: %d\n', num_opt);
fprintf('  Cylinder: V=%.8f m^3, R=%.6f m, H=%.6f m, center=(%.4f, %.4f) m\n', V_opt, R_opt, H_opt, x0_opt, y0_opt);
fprintf('  Max ball-joint angles:  %.2f, %.2f, %.2f, %.2f, %.2f deg\n', max_angles_opt);

if V_opt > V_orig
    improvement = (V_opt - V_orig) / V_orig * 100;
    fprintf('\n>>> Optimized scheme is BETTER, cylinder volume increased by %.2f%% <<<%\n', improvement);
    phi_best = phi_opt;
    angle_data_best = angle_data_opt;
    is_opt_better = true;
else
    fprintf('\n>>> Original scheme is already near-optimal <<<\n');
    phi_best = phi_orig;
    angle_data_best = angle_data_orig;
    is_opt_better = false;
end


% ===== 阶段4：进一步局部精细搜索（以圆柱体积为优化目标） =====
% 在优化方案附近进行局部精细搜索，确认更优解
fprintf('\n>>> Phase 4: Local fine search around optimized phi...\n');

phi_range_deg = 5;   % 局部搜索范围 +/-5°
phi_step_deg = 2.5;  % 步长 2.5°

best_phi_local = phi_opt;
best_V_local = V_opt;

% 逐支链优化策略（避免5维全组合爆炸）
for i_limb = 2 : 5  % limb 1 无角度约束，跳过
    phi_test = best_phi_local;
    phi_center_deg = rad2deg(best_phi_local(i_limb));
    phi_search_deg = max(phi_center_deg - phi_range_deg, 20) : phi_step_deg : min(phi_center_deg + phi_range_deg, 60);
    
    best_phi_for_limb = best_phi_local(i_limb);
    for phi_deg = phi_search_deg
        phi_test(i_limb) = deg2rad(phi_deg);
        [V_test, ~, ~, ~, ~, ~, ~, ~] = eval_scheme(phi_test, theta_vec, seq_x_eval, seq_y_eval, seq_z_eval, P_v, B, l_min, l_max, R_plant, angle_limit_deg);
        if V_test > best_V_local
            best_V_local = V_test;
            best_phi_for_limb = phi_test(i_limb);
        end
    end
    best_phi_local(i_limb) = best_phi_for_limb;
end

if best_V_local > V_opt
    improvement_local = (best_V_local - V_opt) / V_opt * 100;
    fprintf('  Local search found better solution: phi = [%.2f, %.2f, %.2f, %.2f, %.2f] deg\n', rad2deg(best_phi_local));
    fprintf('  Further improvement: %.2f%%\n', improvement_local);
    phi_best = best_phi_local;
    [V_best, R_best, H_best, x0_best, y0_best, num_best, max_angles_best, angle_data_best] = ...
        eval_scheme(phi_best, theta_vec, seq_x_eval, seq_y_eval, seq_z_eval, P_v, B, l_min, l_max, R_plant, angle_limit_deg);
else
    fprintf('  No better solution found in local search.\n');
    if is_opt_better
        phi_best = phi_opt;
        V_best = V_opt;
        R_best = R_opt;
        H_best = H_opt;
        x0_best = x0_opt;
        y0_best = y0_opt;
        num_best = num_opt;
        max_angles_best = max_angles_opt;
    else
        phi_best = phi_orig;
        V_best = V_orig;
        R_best = R_orig;
        H_best = H_orig;
        x0_best = x0_orig;
        y0_best = y0_orig;
        num_best = num_orig;
        max_angles_best = max_angles_orig;
    end
end


% ===== 阶段5：最终方案统计与输出 =====
fprintf('\n========== Final Recommended Scheme ==========\n');
fprintf('Recommended phi (ball_screw_dir_angle row 1):\n');
for i_limb = 1 : 5
    fprintf('  Limb %d: %.4f deg (%.6f rad)\n', i_limb, rad2deg(phi_best(i_limb)), phi_best(i_limb));
end
fprintf('Cylinder volume V = %.8f m^3\n', V_best);
fprintf('Cylinder radius R = %.6f m\n', R_best);
fprintf('Cylinder height H = %.6f m\n', H_best);
fprintf('Cylinder center (x0, y0) = (%.4f, %.4f) m\n', x0_best, y0_best);
fprintf('Diameter/Height ratio = %.6f\n', 2*R_best/H_best);
fprintf('Valid workspace points: %d\n', num_best);
fprintf('Max ball-joint angles:  %.2f, %.2f, %.2f, %.2f, %.2f deg\n', max_angles_best);

% 最终角度统计
fprintf('\nBall-joint angle statistics (final scheme):\n');
for i_limb = 1 : 5
    angles = angle_data_best{i_limb};
    if ~isempty(angles)
        fprintf('  Limb %d: mean=%.2f deg, max=%.2f deg, min=%.2f deg, std=%.2f deg\n', ...
            i_limb, mean(angles), max(angles), min(angles), std(angles));
    end
end

fprintf('\nrotation of plate: x - %.2f deg, y - %.2f deg, z - %.2f deg\n', alpha_plant_deg, beta_plant_deg, gamma_plant_deg);
fprintf('phi (angle to z-axis), theta (angle to x-axis in xoy plane)\n');
for i_limb = 1 : 5
    fprintf('limb %d: phi = %.4f deg, theta = %.4f deg\n', i_limb, rad2deg(phi_best(i_limb)), rad2deg(theta_vec(i_limb)));
end


% ---------------- plot -----------------
if workspace_plot_flag == 1
    % 使用最终方案重新计算工作空间上下边界用于绘图
    seq_x_plot = (-300 : 10 : 300)*unit_para;
    seq_y_plot = (-300 : 10 : 300)*unit_para;
    seq_z_plot = (-1100 : 5 : -700)*unit_para;
    
    ball_vec_best = zeros(3, 5);
    for i = 1 : 5
        ball_vec_best(1, i) = sin(phi_best(i)) * cos(theta_vec(i));
        ball_vec_best(2, i) = sin(phi_best(i)) * sin(theta_vec(i));
        ball_vec_best(3, i) = cos(phi_best(i));
    end
    ball_vec_best = -1 * ball_vec_best;
    ball_vec_world_best = R_plant * ball_vec_best;
    
    work_space_up = [0;0;0];
    work_space_down = [0;0;0];
    
    for ix = 1 : length(seq_x_plot)
        for iy = 1 : length(seq_y_plot)
            z_min_point = zeros(3,1);
            z_max_point = zeros(3,1);
            for iz = 1 : length(seq_z_plot)
                vt = [seq_x_plot(ix); seq_y_plot(iy); seq_z_plot(iz)];
                pos_flag = 0;
                s_limb = zeros(3, 5);
                
                for j = 1 : 5
                    vAa = vt + P_v(:, j) - B(:, j);
                    len_vAa = norm(vAa);
                    s_limb(:, j) = vAa / len_vAa;
                    
                    if (len_vAa >= l_min) && (len_vAa <= l_max)
                        angle_limb_scew = acos(dot(s_limb(:, j), ball_vec_world_best(:, j)));
                        if (rad2deg(angle_limb_scew) <= angle_limit_deg) || (j == 1)
                            pos_flag = pos_flag + 1;
                        end
                    end
                end
                
                if pos_flag == 5
                    if z_min_point == zeros(3,1)
                        z_min_point = vt;
                        z_max_point = vt;
                    else
                        z_max_point = vt;
                    end
                end
            end
            
            if (z_min_point(3) ~= 0) && (z_max_point(3) ~= 0)
                work_space_up = [work_space_up z_max_point];
                work_space_down = [work_space_down z_min_point];
            end
        end
    end
    
    fig = figure('Color', [1 1 1]);
    % 机构简图
    plot3(B(1, :), B(2, :), B(3, :), 'o', 'Color', '#FF7F50');
    hold on
    plot3(P(1, :), P(2, :), P(3, :), 'o', 'Color', '#32CD32');
    B_plot = [B(:,1) B(:,5) B(:,2) B(:,3) B(:,4) B(:,1)];
    P_plot = [P(:,1) P(:,5) P(:,2) P(:,3) P(:,4) P(:,1)];
    plot3(B_plot(1, :), B_plot(2, :), B_plot(3, :), '-', 'Color', '#FF7F50');
    plot3(P_plot(1, :), P_plot(2, :), P_plot(3, :), '-', 'Color', '#32CD32');

    for i = 1 : 5
        plot3([B(1, i) P(1, i)], [B(2, i) P(2, i)], [B(3, i) P(3, i)], '-', 'Color', '#4682B4');
    end
    scatter3(work_space_up(1,2:end), work_space_up(2,2:end), work_space_up(3,2:end), 2, work_space_up(3,2:end),'filled');
    scatter3(work_space_down(1,2:end), work_space_down(2,2:end), work_space_down(3,2:end), 2, work_space_down(3,2:end),'filled');

    grid on
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')

    % view rotation
    axis vis3d
    filename = 'view0.gif';
    fif_delay_time = 0.03;

    if fig_rotation_show == 1
        for ii = 1 : 360
            view(-45 + 1*ii,30);
            pause(fif_delay_time);

            % generate gif
            if gif_generate_flag == 1
                frame = getframe(fig);
                im = frame2im(frame);
                [A, map] = rgb2ind(im, 256);
                if ii == 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',fif_delay_time);
                else
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',fif_delay_time);
                end
            end
        end
    end
end

fprintf('>>>= Joint_angle_search done (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));


%% ============== Local Functions ==============
function [cyl_V, cyl_R, cyl_H, cyl_x0, cyl_y0, num_points, max_angles, limb_angle_data] = ...
    eval_scheme(phi_vec, theta_vec, seq_x, seq_y, seq_z, P_v, B, l_min, l_max, R_plant, angle_limit_deg)
    
    % 构建球铰链方向向量
    ball_vec = zeros(3, 5);
    for i = 1 : 5
        ball_vec(1, i) = sin(phi_vec(i)) * cos(theta_vec(i));
        ball_vec(2, i) = sin(phi_vec(i)) * sin(theta_vec(i));
        ball_vec(3, i) = cos(phi_vec(i));
    end
    ball_vec = -1 * ball_vec;  % 与主脚本一致，取反使方向与支链同向
    ball_vec_world = R_plant * ball_vec;
    
    len_x = length(seq_x);
    len_y = length(seq_y);
    len_z = length(seq_z);
    is_reachable = false(len_x, len_y, len_z);
    
    num_points = 0;
    max_angles = zeros(1, 5);
    limb_angle_data = cell(1, 5);
    for i = 1 : 5
        limb_angle_data{i} = [];
    end
    
    % ===== 遍历网格，搜索可达点并构建 is_reachable =====
    for ix = 1 : len_x
        for iy = 1 : len_y
            for iz = 1 : len_z
                vt = [seq_x(ix); seq_y(iy); seq_z(iz)];
                pos_flag = 0;
                s_limb = zeros(3, 5);
                limb_angles = zeros(1, 5);
                
                for j = 1 : 5
                    vAa = vt + P_v(:, j) - B(:, j);
                    len_vAa = norm(vAa);
                    s_limb(:, j) = vAa / len_vAa;
                    
                    if (len_vAa >= l_min) && (len_vAa <= l_max)
                        angle_limb_scew = acos(dot(s_limb(:, j), ball_vec_world(:, j)));
                        limb_angles(j) = rad2deg(angle_limb_scew);
                        % 与v3一致：limb 1 免除角度约束，其余 limb 需满足 <= 30°
                        if (limb_angles(j) <= angle_limit_deg) || (j == 1)
                            pos_flag = pos_flag + 1;
                        end
                    end
                end
                
                if pos_flag == 5
                    is_reachable(ix, iy, iz) = true;
                    num_points = num_points + 1;
                    max_angles = max(max_angles, limb_angles);
                    for j = 1 : 5
                        limb_angle_data{j} = [limb_angle_data{j} limb_angles(j)];
                    end
                end
            end
        end
    end
    
    % ===== 圆柱检索（参照 workspace_discrete_v3.m） =====
    % 约束参数
    H_limit = 0.08;     % 最小高度 0.08 m
    R_limit = 0.055;    % 最小半径 0.055 m
    ratio_limit = 0.6;  % 最小直径/高度比 2R/H > 0.6
    
    cyl_V = 0;
    cyl_R = 0;
    cyl_H = 0;
    cyl_x0 = 0;
    cyl_y0 = 0;
    
    for iz1 = 1 : len_z
        for iz2 = iz1 : len_z
            H_cand = seq_z(iz2) - seq_z(iz1);
            if H_cand < H_limit
                continue;
            end
            
            % 该 z 区间内的共同可达点
            common_ws = squeeze(is_reachable(:, :, iz1));
            for iz = iz1+1 : iz2
                common_ws = common_ws & squeeze(is_reachable(:, :, iz));
            end
            
            if ~any(common_ws(:))
                continue;
            end
            
            % 可达点与不可达点坐标
            [reach_ix, reach_iy] = find(common_ws);
            reach_pts = [seq_x(reach_ix).', seq_y(reach_iy).'];
            
            [unreach_ix, unreach_iy] = find(~common_ws);
            unreach_pts = [seq_x(unreach_ix).', seq_y(unreach_iy).'];
            
            % 计算该 z 区间内的最大内接圆（圆心任意）
            best_R_for_interval = 0;
            best_center_for_interval = [0, 0];
            
            if isempty(unreach_pts)
                % 全部可达，圆心取质心，半径为最远的可达点距离
                x0_c = mean(reach_pts(:,1));
                y0_c = mean(reach_pts(:,2));
                dists = sqrt((reach_pts(:,1)-x0_c).^2 + (reach_pts(:,2)-y0_c).^2);
                best_R_for_interval = max(dists);
                best_center_for_interval = [x0_c, y0_c];
            else
                % 使用 bsxfun 批量计算可达点到不可达点的距离
                dx = bsxfun(@minus, reach_pts(:,1), unreach_pts(:,1)');
                dy = bsxfun(@minus, reach_pts(:,2), unreach_pts(:,2)');
                D = sqrt(dx.^2 + dy.^2);
                R_limits = min(D, [], 2);  % 每个候选圆心的最近不可达点距离
                
                for i = 1 : size(reach_pts, 1)
                    R_safe = R_limits(i) - 1e-6;
                    if R_safe > best_R_for_interval
                        best_R_for_interval = R_safe;
                        best_center_for_interval = reach_pts(i, :);
                    end
                end
            end
            
            R = best_R_for_interval;
            if R < R_limit
                continue;
            end
            if (2*R/H_cand) <= ratio_limit
                continue;
            end
            
            V = pi * R^2 * H_cand;
            if V > cyl_V
                cyl_V = V;
                cyl_R = R;
                cyl_H = H_cand;
                cyl_x0 = best_center_for_interval(1);
                cyl_y0 = best_center_for_interval(2);
            end
        end
    end
end
