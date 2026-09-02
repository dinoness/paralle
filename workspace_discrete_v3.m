% 计算满足一定角度的圆柱工作空间，并计算该圆柱空间的内的OTI，LCI
% >>>>>所描述的空间是刀具末端点<<<<<
% =====================注意搜索间隙=====================
clear
fig_rotation_show = 0;  % 1开启展示旋转
gif_generate_flag = 0;  % 1为开启录制功能，运行一次程序后记得改文件名
flag_range_plot = 0;  % 1为开启，绘制出指定高度的空间范围

% para_control
H_limit_red = 0.08;
R_limit_red = 0.055;
H_limit_blue = 0.1;
R_limit_blue = 0.1;

path_add()
fprintf('>>>= start (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));
%--------parameter3--------
basic_paras = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'm');  % 单位意思是程序中用到的单位

% >>> 自动加载 surrogateopt 优化结果 <<<
use_optimized_params = false;

if use_optimized_params && exist('optimization_result.mat', 'file')
    load('optimization_result.mat', 'results');
    x_opt = results.x_opt(:)';
    fprintf('\n>>> Loaded optimized parameters from optimization_result.mat <<<\n');
    fprintf('    x_opt = [H=%.4f mm, h=%.4f mm, r1=%.4f mm, r2=%.4f mm, a=%.4f deg, b=%.4f deg]\n', ...
            x_opt(1), x_opt(2), x_opt(3), x_opt(4), x_opt(5), x_opt(6));

    % 按 evaluate_spr4ups_objective.m 的相同方式覆盖参数
    unit_para = basic_paras.unit_para;
    l_max = basic_paras.l_max;
    l_min = basic_paras.l_min;
    R1    = basic_paras.R1;
    R2    = basic_paras.R2;
    H     = x_opt(1) * 0.001;
    h     = x_opt(2) * 0.001;
    r1    = x_opt(3) * 0.001;
    r2    = x_opt(4) * 0.001;
    L_tool = basic_paras.L_tool;
    joint_u_angle_tilt = basic_paras.joint_u_angle_tilt;
    l0_seq = basic_paras.l0_seq;

    % 更新 limb_dir（与 evaluate_spr4ups_objective.m 一致）
    limb_dir = basic_paras.limb_dir;
    a = deg2rad(x_opt(5));
    b = deg2rad(x_opt(6));
    limb_dir(4, 2) = deg2rad(90) - a;
    limb_dir(5, 2) = deg2rad(90) + a;
    limb_dir(2, 2) = deg2rad(270) - b;
    limb_dir(3, 2) = deg2rad(270) + b;

    % 重新计算 B, P_m（与 evaluate_spr4ups_objective.m 一致）
    B1 = [R1*cos(limb_dir(1,1)); R1*sin(limb_dir(1,1)); 0];
    B2 = [R1*cos(limb_dir(2,1)); R1*sin(limb_dir(2,1)); 0];
    B3 = [R1*cos(limb_dir(3,1)); R1*sin(limb_dir(3,1)); 0];
    B4 = [R2*cos(limb_dir(4,1)); R2*sin(limb_dir(4,1)); H];
    B5 = [R2*cos(limb_dir(5,1)); R2*sin(limb_dir(5,1)); H];
    B = [B1 B2 B3 B4 B5];

    P1_m = [r1*cos(limb_dir(1,2)); r1*sin(limb_dir(1,2)); L_tool];
    P2_m = [r1*cos(limb_dir(2,2)); r1*sin(limb_dir(2,2)); L_tool];
    P3_m = [r1*cos(limb_dir(3,2)); r1*sin(limb_dir(3,2)); L_tool];
    P4_m = [r2*cos(limb_dir(4,2)); r2*sin(limb_dir(4,2)); L_tool+h];
    P5_m = [r2*cos(limb_dir(5,2)); r2*sin(limb_dir(5,2)); L_tool+h];
    P_m = [P1_m P2_m P3_m P4_m P5_m];

    fprintf('>>> Optimized parameters applied successfully. <<<\n\n');
else
    % 使用 parameters.xlsx 原始参数
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

    if use_optimized_params
        warning('optimization_result.mat not found. Using parameters.xlsx directly.');
    end
end

pos_plant = [0; 0; -800]*unit_para;  % 后面作图用，不参与空间搜索
R_plant = eye(3);
P_v = zeros(3, 5);  % 只变换了方向，没变换起点
P = zeros(3, 5);    % 末端点坐标
for i = 1 : 5
    P_v(:, i) = R_plant * P_m(:, i);
    P(:, i) = P_v(:, i) + pos_plant;
end


% ball screw dir vector
ball_screw_dir_angle = [deg2rad([35 35 35 35 35]);
                        limb_dir(:, 2)'];
ball_vector = zeros(3, 5);
ball_vector_world = zeros(3, 5);

% 目前未考虑虎克铰
static_joint_dir_angle_deg = [ 0  0  0  0  0;
                              -90  30 150 -135 -45];  
static_joint_dir_angle = static_joint_dir_angle_deg / 180 * pi;
static_joint_vector = zeros(3, 5);

for i  = 1 : 5
    ball_vector(1, i ) = sin(ball_screw_dir_angle(1, i )) * cos(ball_screw_dir_angle(2, i ));
    ball_vector(2, i ) = sin(ball_screw_dir_angle(1, i )) * sin(ball_screw_dir_angle(2, i ));
    ball_vector(3, i ) = cos(ball_screw_dir_angle(1, i ));

    static_joint_vector(1, i ) = sin(static_joint_dir_angle(1, i )) * cos(static_joint_dir_angle(2, i ));
    static_joint_vector(2, i ) = sin(static_joint_dir_angle(1, i )) * sin(static_joint_dir_angle(2, i ));
    static_joint_vector(3, i ) = cos(static_joint_dir_angle(1, i ));
end
ball_vector = -1 * ball_vector;  % 让关节方向向量和支链方向同向

% 预计算 parameterize 参数，供 LTI/GCI 使用
p_seq = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

% -----end-parameter3------



% ------search space-------
seq_x = (-400 : 10 : 400)*unit_para;
seq_y = (-400 : 10 : 400)*unit_para;
seq_z = (-1200 : 5 : -600)*unit_para;
seq_phi = deg2rad((-180 : 30 : 180));
seq_theta = deg2rad((0 : 2 : 20));
len_x = length(seq_x);
len_y = length(seq_y);
len_z = length(seq_z);
len_theta = length(seq_theta);
len_phi = length(seq_phi);


% assistant parameter
work_space_ang = [0;0;0;0];  % 含角度的工作空间
ang_threshold = deg2rad(9);
work_space_up = [0;0;0];
work_space_down = [0;0;0];
pos_count = 0;  % 空间点计数

% 记录每个 [x,y,z] 下哪些 theta 是可达的（所有 phi 均满足）
reachable_thetas = false(len_x, len_y, len_z, len_theta);
pos_quality_matrix = -ones(len_x, len_y, len_z);  % 记录每个点的 pos_quality（可达的最大theta），-1为不可达


%% 遍历
parfor ix = 1 : len_x
    px = seq_x(ix);
    local_ws_ang = [0;0;0;0];
    local_ws_up = [0;0;0];
    local_ws_down = [0;0;0];
    local_reachable = false(len_y, len_z, len_theta);
    local_pos_quality = -ones(len_y, len_z);
    local_count = 0;

    for iy = 1 : len_y
        py = seq_y(iy);
        % x,y方向进行遍历

        z_min_point = zeros(3,1);
        z_max_point = zeros(3,1);

        % 搜索z向上的工作空间界限
        for iz = 1 : len_z
            pz = seq_z(iz);
            pos_quality = -1;  % 代表改点下姿态可达性，-1为改点不可达，值表示其可达的theta角度

            for itheta = 1 : len_theta
                flag_all_phi = 0;  % 记录在指定theta下，是否phi都可达
                theta = seq_theta(itheta);

                for iphi = 1 : len_phi
                    pos_flag = 0;  % 位置可达标志位
                    s_limb = zeros(3, 5);  % 支链的方向向量
                    l_limb = zeros(1, 5);  % 支链长度

                    Pos_ref = [px; py; pz; seq_phi(iphi); theta];
                    T_ref = pos2trans(Pos_ref, B, 'unit', 'rad');
                    vt = T_ref(1:3, 4);
                    R_plant = T_ref(1:3, 1:3);
                    P_v = R_plant * P_m;
                    ball_vector_world = R_plant * ball_vector;
                    

                    for j = 1 : length(P_v(1, :))
                        vAa = vt + P_v(:, j) - B(:, j);
                        len_vAa = norm(vAa);
                        l_limb(j) = len_vAa;
                        s_limb(:, j) = vAa / len_vAa;

                        if (len_vAa >= l_min)&&(len_vAa <= l_max)  % ===========支链长度条件===========
                            angle_limb_scew = acos(dot(s_limb(:, j), ball_vector_world(:, j)));  % 支链与球铰轴线夹角
                            angle_limb_scew_deg = rad2deg(angle_limb_scew);
                            if(angle_limb_scew_deg <= 30 ) || (j == 1)  % ===========关节角度条件===========
                                pos_flag = pos_flag + 1;
                            end
                        end
                    end

                    if(pos_flag < length(P_v(1, :)))
                        flag_all_phi = -1;
                        break;
                    end

                end  % phi

                % 记录theta情况
                if(flag_all_phi < 0)
                    break;
                else
                    pos_quality = theta;
                    local_reachable(iy, iz, itheta) = true;
                end
            end  % theta

            local_pos_quality(iy, iz) = pos_quality;

            % 含角度的工作空间
            if pos_quality > ang_threshold
                p_mark = [vt; rad2deg(pos_quality)];
                local_ws_ang = [local_ws_ang p_mark];
            end

            % 工作空间轮廓：只画z方向的上下两端点
            if pos_quality > -1  % 如果所有条件均允许
                if z_min_point == zeros(3,1)  % 如果第一次进入循环
                    z_min_point = vt;
                    z_max_point = vt;
                else
                    z_max_point = vt;
                end

                local_count = local_count + 1;
            end

        end  % z

        % 工作空间轮廓：将工作空间界限添加到最后的作图中
        if (z_min_point(3) ~= 0) && (z_max_point(3) ~= 0)
            local_ws_up = [local_ws_up z_max_point];
            local_ws_down = [local_ws_down z_min_point];
        end

    end  % y

    % 合并到全局变量（parfor reduction）
    work_space_ang = [work_space_ang local_ws_ang(:, 2:end)];
    work_space_up = [work_space_up local_ws_up(:, 2:end)];
    work_space_down = [work_space_down local_ws_down(:, 2:end)];
    pos_count = pos_count + local_count;
    reachable_thetas(ix, :, :, :) = local_reachable;
    pos_quality_matrix(ix, :, :) = local_pos_quality;
end  % x

fprintf('Workspace search done. Reachable points with theta > threshold: %d\n', size(work_space_ang, 2)-1);



%% 含角度的工作空间
fig_ang = figure('Color', [1 1 1]);
scatter3(work_space_ang(1,2:end), work_space_ang(2,2:end), work_space_ang(3,2:end), 2, work_space_ang(4,2:end), 'filled');
grid on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
colormap(parula);
colorbar;


%% ================== Cylinder Extraction ==================
% 从 work_space_ang 中提取三维逻辑矩阵
% work_space_ang 的每一列为 [x; y; z; theta_max_deg]
is_ws_ang = false(len_x, len_y, len_z);
for k = 2 : size(work_space_ang, 2)
    [~, ix] = min(abs(seq_x - work_space_ang(1, k)));
    [~, iy] = min(abs(seq_y - work_space_ang(2, k)));
    [~, iz] = min(abs(seq_z - work_space_ang(3, k)));
    is_ws_ang(ix, iy, iz) = true;
end

% 寻找满足约束的最大体积圆柱
% 约束：R >= 0.05 m, H >= 0.15 m, 2R/H > 0.6
% 圆柱底面平行 xOy，轴线平行 z，圆心 (x0, y0) 任意
% 圆柱必须完全包含在 work_space_ang 的离散点集中
best_V = 0;
best_R = 0;
best_H = 0;
best_x0 = 0;
best_y0 = 0;
best_iz1 = 0;
best_iz2 = 0;

for iz1 = 1 : len_z
    for iz2 = iz1 : len_z
        H = seq_z(iz2) - seq_z(iz1);
        if H < H_limit_red
            continue;
        end
        
        % 该 z 区间内的共同可达点（在所有 z 层上都是 work_space_ang）
        common_ws = squeeze(is_ws_ang(:, :, iz1));
        for iz = iz1+1 : iz2
            common_ws = common_ws & squeeze(is_ws_ang(:, :, iz));
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
        if R < R_limit_red
            continue;
        end
        if (2*R/H) <= 0.6
            continue;
        end
        
        V = pi * R^2 * H;
        if V > best_V
            best_V = V;
            best_R = R;
            best_H = H;
            best_x0 = best_center_for_interval(1);
            best_y0 = best_center_for_interval(2);
            best_iz1 = iz1;
            best_iz2 = iz2;
        end
    end
end

if best_V <= 0
    warning('No cylinder satisfies the constraints inside work_space_ang (R>=%.4f, H>=%.4f, 2R/H>0.6).', R_limit_red,H_limit_red);
    fprintf('Diagnostics:\n');
    fprintf('  work_space_ang points (theta > %.1f deg): %d\n', rad2deg(ang_threshold), size(work_space_ang,2)-1);
    if size(work_space_ang,2) > 1
        z_unique = unique(work_space_ang(3,2:end));
        fprintf('  Unique z layers in work_space_ang: %d (range: [%.4f, %.4f] m)\n', ...
                length(z_unique), min(z_unique), max(z_unique));
    end
    fprintf('Skipping LTI/GCI and LCI computations for cylinder 1.\n\n');

    R_cyl = 0; H_cyl = 0; x0_cyl = 0; y0_cyl = 0;
    z_min_cyl = 0; z_max_cyl = 0; V_cyl = 0;
    best_iz1 = 1; best_iz2 = 1;
else
    R_cyl = best_R;
    H_cyl = best_H;
    x0_cyl = best_x0;
    y0_cyl = best_y0;
    z_min_cyl = seq_z(best_iz1);
    z_max_cyl = seq_z(best_iz2);
    V_cyl = best_V;

    fprintf('\n========== Cylinder Workspace ==========\n');
    fprintf('Cylinder center (x0, y0) = (%.6f, %.6f) m\n', x0_cyl, y0_cyl);
    fprintf('Cylinder radius R = %.6f m\n', R_cyl);
    fprintf('Cylinder height H = %.6f m\n', H_cyl);
    fprintf('Diameter/Height ratio = %.6f\n', 2*R_cyl/H_cyl);
    fprintf('Cylinder volume V = %.8f m^3\n', V_cyl);
    fprintf('z range = [%.4f, %.4f] m\n', z_min_cyl, z_max_cyl);
    fprintf('========================================\n\n');
end


if best_V > 0
    %% ================== LTI / GCI Computation ==================
    fprintf('Starting LTI/GCI computation... This may take a while.\n');

    % 分块并行计算，减少通信开销
    ltigci_cell = cell(len_x, 1);

    parfor ix = 1 : len_x
        px = seq_x(ix);
        local_data = [];
        for iy = 1 : len_y
            py = seq_y(iy);
            r_xy = sqrt((px - x0_cyl)^2 + (py - y0_cyl)^2);
            if r_xy > R_cyl
                continue;
            end
            for iz = 1 : len_z
                pz = seq_z(iz);
                if pz < z_min_cyl || pz > z_max_cyl
                    continue;
                end
                for itheta = 1 : len_theta
                    if ~reachable_thetas(ix, iy, iz, itheta)
                        continue;
                    end
                    theta = seq_theta(itheta);

                    lti_phi = zeros(len_phi, 1);
                    gci_phi = zeros(len_phi, 1);
                    for iphi = 1 : len_phi
                        Pos_ref = [px; py; pz; seq_phi(iphi); theta];
                        T_ref = pos2trans(Pos_ref, B, 'unit', 'rad');
                        [lti_val, gci_val] = compute_ltigci(T_ref, B, l0_seq, P_m, p_seq);
                        lti_phi(iphi) = lti_val;
                        gci_phi(iphi) = gci_val;
                    end
                    mean_lti = mean(lti_phi);
                    mean_gci = mean(gci_phi);
                    local_data = [local_data; px, py, pz, theta, mean_lti, mean_gci, ix, iy, iz];
                end
            end
        end
        ltigci_cell{ix} = local_data;
    end

    % 合并结果
    all_ltigci = [];
    for ix = 1 : len_x
        if ~isempty(ltigci_cell{ix})
            all_ltigci = [all_ltigci; ltigci_cell{ix}];
        end
    end

    N_total = size(all_ltigci, 1);
    if N_total == 0
        warning('No reachable points found inside the cylinder. Skipping LTI/GCI output.');
    else
        % 输出 1：全局统计
        mean_LTI = mean(all_ltigci(:,5));
        min_LTI = min(all_ltigci(:,5));
        mean_GCI = mean(all_ltigci(:,6));
        min_GCI = min(all_ltigci(:,6));

        fprintf('\n========== LTI / GCI Statistics (all [x,y,z,theta] in cylinder) ==========\n');
        fprintf('Total evaluated [x,y,z,theta] points: %d\n', N_total);
        fprintf('LTI  —  mean: %.8f,  min: %.8f\n', mean_LTI, min_LTI);
        fprintf('GCI  —  mean: %.8f,  min: %.8f\n', mean_GCI, min_GCI);
        fprintf('==========================================================================\n\n');

        %% ================== Output 2: 3D Scatter Plot ==================
        % 对每个 [x,y,z]，在可达的 theta 上取 LTI 和 GCI 的平均值
        idx_x = all_ltigci(:, 7);
        idx_y = all_ltigci(:, 8);
        idx_z = all_ltigci(:, 9);

        lin_idx = sub2ind([len_x, len_y, len_z], idx_x, idx_y, idx_z);
        lti_sum = accumarray(lin_idx, all_ltigci(:,5), [len_x*len_y*len_z, 1]);
        lti_cnt = accumarray(lin_idx, ones(N_total,1), [len_x*len_y*len_z, 1]);
        gci_sum = accumarray(lin_idx, all_ltigci(:,6), [len_x*len_y*len_z, 1]);
        gci_cnt = accumarray(lin_idx, ones(N_total,1), [len_x*len_y*len_z, 1]);

        valid_idx = find(lti_cnt > 0);
        N_valid = length(valid_idx);
        point_LTI = zeros(N_valid, 1);
        point_GCI = zeros(N_valid, 1);
        point_x = zeros(N_valid, 1);
        point_y = zeros(N_valid, 1);
        point_z = zeros(N_valid, 1);

        for k = 1 : N_valid
            lin = valid_idx(k);
            [ix, iy, iz] = ind2sub([len_x, len_y, len_z], lin);
            point_x(k) = seq_x(ix);
            point_y(k) = seq_y(iy);
            point_z(k) = seq_z(iz);
            point_LTI(k) = lti_sum(lin) / lti_cnt(lin);
            point_GCI(k) = gci_sum(lin) / gci_cnt(lin);
        end

        fig_cyl = figure('Color', [1 1 1]);
        scatter3(point_x, point_y, point_z, 20, point_LTI, 'filled');
        % scatter3(work_space_ang(1,2:end), work_space_ang(2,2:end), work_space_ang(3,2:end), 2, work_space_ang(4,2:end), 'filled');
        grid on
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')
        colormap(jet);
        colorbar;
        title('Cylindrical Workspace colored by LTI (averaged over \theta)');

        % 绘制圆柱边界（可选，帮助可视化）
        hold on
        % 画上下两个圆
        theta_circle = linspace(0, 2*pi, 100);
        x_circle_top = x0_cyl + R_cyl * cos(theta_circle);
        y_circle_top = y0_cyl + R_cyl * sin(theta_circle);
        z_circle_top = z_max_cyl * ones(size(theta_circle));
        z_circle_bottom = z_min_cyl * ones(size(theta_circle));
        plot3(x_circle_top, y_circle_top, z_circle_top, 'r--', 'LineWidth', 1.5);
        plot3(x_circle_top, y_circle_top, z_circle_bottom, 'r--', 'LineWidth', 1.5);
        for k_line = 1 : 4 : length(theta_circle)
            plot3([x_circle_top(k_line), x_circle_top(k_line)], ...
                  [y_circle_top(k_line), y_circle_top(k_line)], ...
                  [z_min_cyl, z_max_cyl], 'r--', 'LineWidth', 0.5);
        end
        hold off
    end
end



%% ================== 可达空间中的最大圆柱 (R>0.1, H>0.1) ==================
% 基于所有可达点（pos_quality >= 0）
is_reachable_anytheta = pos_quality_matrix >= 0;

% 寻找满足 R > 0.1, H > 0.1 的最大体积圆柱
best_V2 = 0;
best_R2 = 0;
best_H2 = 0;
best_x02 = 0;
best_y02 = 0;
best_iz1_2 = 0;
best_iz2_2 = 0;

for iz1 = 1 : len_z
    for iz2 = iz1 : len_z
        H_cand = seq_z(iz2) - seq_z(iz1);
        if H_cand <= H_limit_blue
            continue;
        end
        
        common_r = squeeze(is_reachable_anytheta(:, :, iz1));
        for iz = iz1+1 : iz2
            common_r = common_r & squeeze(is_reachable_anytheta(:, :, iz));
        end
        
        if ~any(common_r(:))
            continue;
        end
        
        [reach_ix, reach_iy] = find(common_r);
        reach_pts = [seq_x(reach_ix).', seq_y(reach_iy).'];
        
        [unreach_ix, unreach_iy] = find(~common_r);
        unreach_pts = [seq_x(unreach_ix).', seq_y(unreach_iy).'];
        
        best_R_for_interval = 0;
        best_center_for_interval = [0, 0];
        
        if isempty(unreach_pts)
            x0_c = mean(reach_pts(:,1));
            y0_c = mean(reach_pts(:,2));
            dists = sqrt((reach_pts(:,1)-x0_c).^2 + (reach_pts(:,2)-y0_c).^2);
            best_R_for_interval = max(dists);
            best_center_for_interval = [x0_c, y0_c];
        else
            dx = bsxfun(@minus, reach_pts(:,1), unreach_pts(:,1)');
            dy = bsxfun(@minus, reach_pts(:,2), unreach_pts(:,2)');
            D = sqrt(dx.^2 + dy.^2);
            R_limits = min(D, [], 2);
            
            for i = 1 : size(reach_pts, 1)
                R_safe = R_limits(i) - 1e-6;
                if R_safe > best_R_for_interval
                    best_R_for_interval = R_safe;
                    best_center_for_interval = reach_pts(i, :);
                end
            end
        end
        
        R = best_R_for_interval;
        if R <= R_limit_blue
            continue;
        end
        if (2*R/H_cand) <= 0.6
            continue;
        end
        
        V = pi * R^2 * H_cand;
        if V > best_V2
            best_V2 = V;
            best_R2 = R;
            best_H2 = H_cand;
            best_x02 = best_center_for_interval(1);
            best_y02 = best_center_for_interval(2);
            best_iz1_2 = iz1;
            best_iz2_2 = iz2;
        end
    end
end

fprintf('\n========== Cylinder in Reachable Space (R>0.1, H>0.1) ==========\n');
if best_V2 <= 0
    fprintf('No cylinder satisfies R > %.4f and H > %.4f in reachable space.\n', R_limit_blue, H_limit_blue);
else
    fprintf('Center (x0, y0) = (%.6f, %.6f) m\n', best_x02, best_y02);
    fprintf('Radius R = %.6f m\n', best_R2);
    fprintf('Height H = %.6f m\n', best_H2);
    fprintf('Diameter/Height ratio = %.6f\n', 2*best_R2/best_H2);
    fprintf('Volume V = %.8f m^3\n', best_V2);
    fprintf('z range = [%.4f, %.4f] m\n', seq_z(best_iz1_2), seq_z(best_iz2_2));
end
fprintf('================================================================\n\n');

%% ================== OTI at phi=0, theta=0 in the cylinder ==================
if best_V2 > 0
    oti_vals = [];
    for ix = 1 : len_x
        for iy = 1 : len_y
            px = seq_x(ix);
            py = seq_y(iy);
            if sqrt((px-best_x02)^2 + (py-best_y02)^2) > best_R2 + 1e-9
                continue;
            end
            for iz = best_iz1_2 : best_iz2_2
                if ~is_reachable_anytheta(ix, iy, iz)
                    continue;
                end
                pz = seq_z(iz);
                % 姿态 phi=0, theta=0
                Pos_ref = [px; py; pz; 0; 0];
                T_ref = pos2trans(Pos_ref, B, 'unit', 'rad');
                [~, ~, oti_val, ~] = compute_ltigci(T_ref, B, l0_seq, P_m, p_seq);
                oti_vals = [oti_vals; oti_val];
            end
        end
    end
    
    if ~isempty(oti_vals)
        fprintf('========== OTI at phi=0, theta=0 in the cylinder ==========\n');
        fprintf('Evaluated points: %d\n', length(oti_vals));
        fprintf('OTI mean = %.8f\n', mean(oti_vals));
        fprintf('OTI min  = %.8f\n', min(oti_vals));
        fprintf('OTI max  = %.8f\n', max(oti_vals));
        fprintf('============================================================\n\n');
    else
        fprintf('No valid OTI points at phi=0, theta=0 inside the cylinder.\n');
    end
end

%% ================== LCI 计算（两个圆柱空间） ==================
% LCI = 1 / cond(J)，其中 J 为雅可比矩阵
% 参考 workspace_discrete_v2.m 中的条件数计算方法

% 圆柱 1：基于 work_space_ang 的圆柱
% 遍历圆柱内所有可达的 [x,y,z,theta,phi] 计算 LCI
lci_cyl1 = [];
if best_V > 0
    for ix = 1 : len_x
        for iy = 1 : len_y
            px = seq_x(ix);
            py = seq_y(iy);
            if sqrt((px-x0_cyl)^2 + (py-y0_cyl)^2) > R_cyl + 1e-9
                continue;
            end
            for iz = best_iz1 : best_iz2
                pz = seq_z(iz);
                if ~is_ws_ang(ix, iy, iz)
                    continue;
                end
                for itheta = 1 : len_theta
                    if ~reachable_thetas(ix, iy, iz, itheta)
                        continue;
                    end
                    theta = seq_theta(itheta);
                    for iphi = 1 : len_phi
                        phi = seq_phi(iphi);
                        Pos_ref = [px; py; pz; phi; theta];
                        T_ref = pos2trans(Pos_ref, B, 'unit', 'rad');
                        vt = T_ref(1:3, 4);
                        R_plant = T_ref(1:3, 1:3);
                        P_v = R_plant * P_m;

                        s_limb = zeros(3, 5);
                        l_limb = zeros(1, 5);
                        for j = 1 : 5
                            vAa = vt + P_v(:, j) - B(:, j);
                            len_vAa = norm(vAa);
                            l_limb(j) = len_vAa;
                            s_limb(:, j) = vAa / len_vAa;
                        end

                        x_m = [1; 0; 0];
                        J1 = [s_limb (R_plant * x_m)];
                        J2 = [cross(P_v, s_limb, 1)  ...
                              (cross(P_v(:, 1), (R_plant * x_m)) + l_limb(1)*cross((R_plant * x_m), s_limb(:, 1)))];
                        J = [J1' J2'];
                        lci_cyl1 = [lci_cyl1; 1/cond(J)];
                    end
                end
            end
        end
    end
end

% 圆柱 2：基于所有可达点的圆柱
lci_cyl2 = [];
if best_V2 > 0
    for ix = 1 : len_x
        for iy = 1 : len_y
            px = seq_x(ix);
            py = seq_y(iy);
            if sqrt((px-best_x02)^2 + (py-best_y02)^2) > best_R2 + 1e-9
                continue;
            end
            for iz = best_iz1_2 : best_iz2_2
                pz = seq_z(iz);
                if ~is_reachable_anytheta(ix, iy, iz)
                    continue;
                end
                Pos_ref = [px; py; pz; 0; 0];
                T_ref = pos2trans(Pos_ref, B, 'unit', 'rad');
                vt = T_ref(1:3, 4);
                R_plant = T_ref(1:3, 1:3);
                P_v = R_plant * P_m;

                s_limb = zeros(3, 5);
                l_limb = zeros(1, 5);
                for j = 1 : 5
                    vAa = vt + P_v(:, j) - B(:, j);
                    len_vAa = norm(vAa);
                    l_limb(j) = len_vAa;
                    s_limb(:, j) = vAa / len_vAa;
                end

                x_m = [1; 0; 0];
                J1 = [s_limb (R_plant * x_m)];
                J2 = [cross(P_v, s_limb, 1)  ...
                      (cross(P_v(:, 1), (R_plant * x_m)) + l_limb(1)*cross((R_plant * x_m), s_limb(:, 1)))];
                J = [J1' J2'];
                lci_cyl2 = [lci_cyl2; 1/cond(J)];
            end
        end
    end
end

fprintf('\n========== LCI Statistics ==========\n');
if ~isempty(lci_cyl1)
    fprintf('Cylinder 1 (work_space_ang, all theta&phi): mean LCI = %.8f, points = %d\n', mean(lci_cyl1), length(lci_cyl1));
else
    fprintf('Cylinder 1: no valid points for LCI.\n');
end
if ~isempty(lci_cyl2)
    fprintf('Cylinder 2 (reachable, phi=0, theta=0):     mean LCI = %.8f, points = %d\n', mean(lci_cyl2), length(lci_cyl2));
else
    fprintf('Cylinder 2: no valid points for LCI.\n');
end
fprintf('====================================\n\n');


%% ------- plot workspace (same as v2) -------
fig = figure('Color', [1 1 1]);
% 机构简图
plot3(B(1, :), B(2, :), B(3, :), 'o', 'Color', '#FF7F50', 'LineWidth', 1.5);
hold on
plot3(P(1, :), P(2, :), P(3, :), 'o', 'Color', '#32CD32', 'LineWidth', 1.5);
B_plot = [B(:,1) B(:,5) B(:,2) B(:,3) B(:,4) B(:,1)];
P_plot = [P(:,1) P(:,5) P(:,2) P(:,3) P(:,4) P(:,1)];
plot3(B_plot(1, :), B_plot(2, :), B_plot(3, :), '-', 'Color', '#FF7F50', 'LineWidth', 1.5);
plot3(P_plot(1, :), P_plot(2, :), P_plot(3, :), '-', 'Color', '#32CD32', 'LineWidth', 1.5);

for i = 1 : 5
    plot3([B(1, i) P(1, i)], [B(2, i) P(2, i)], [B(3, i) P(3, i)], '-', 'Color', '#4682B4', 'LineWidth', 1.5);
end
scatter3(work_space_up(1,2:end), work_space_up(2,2:end), work_space_up(3,2:end), 2, work_space_up(3,2:end),'filled');
scatter3(work_space_down(1,2:end), work_space_down(2,2:end), work_space_down(3,2:end), 2, work_space_down(3,2:end),'filled');

if best_V > 0
    plot3(x_circle_top, y_circle_top, z_circle_top, 'r--', 'LineWidth', 1.5);
    plot3(x_circle_top, y_circle_top, z_circle_bottom, 'r--', 'LineWidth', 1.5);
    for k_line = 1 : 4 : length(theta_circle)
        plot3([x_circle_top(k_line), x_circle_top(k_line)], ...
              [y_circle_top(k_line), y_circle_top(k_line)], ...
              [z_min_cyl, z_max_cyl], 'r--', 'LineWidth', 0.5);
    end
end

% 绘制新圆柱（可达空间中的最大圆柱 R>0.1, H>0.1），用蓝色虚线区分
if best_V2 > 0
    theta_circle2 = linspace(0, 2*pi, 100);
    x_circle2 = best_x02 + best_R2 * cos(theta_circle2);
    y_circle2 = best_y02 + best_R2 * sin(theta_circle2);
    z_top2 = seq_z(best_iz2_2) * ones(size(theta_circle2));
    z_bot2 = seq_z(best_iz1_2) * ones(size(theta_circle2));
    plot3(x_circle2, y_circle2, z_top2, 'b--', 'LineWidth', 1.5);
    plot3(x_circle2, y_circle2, z_bot2, 'b--', 'LineWidth', 1.5);
    for k_line = 1 : 4 : length(theta_circle2)
        plot3([x_circle2(k_line), x_circle2(k_line)], ...
              [y_circle2(k_line), y_circle2(k_line)], ...
              [seq_z(best_iz1_2), seq_z(best_iz2_2)], 'b--', 'LineWidth', 0.5);
    end
end
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'LineWidth', 1.5);

legend_str = {};
if best_V > 0
    legend_str{end+1} = 'Workspace cylinder (work_space_ang, R≥0.055)';
end
if best_V2 > 0
    legend_str{end+1} = 'Workspace cylinder (reachable, R>0.1)';
end
if ~isempty(legend_str)
    legend(legend_str, 'Location', 'best');
end
hold off

grid on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

fprintf('>>>= done (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));
