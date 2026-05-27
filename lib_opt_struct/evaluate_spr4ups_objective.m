function J = evaluate_spr4ups_objective(params, k1, k2)
%EVALUATE_SPR4UPS_OBJECTIVE 计算 SPR-4UPS 结构参数下的优化目标值
%
%   输入：
%     params — 6×1 向量 [H_mm, h_mm, r1_mm, r2_mm, a_deg, b_deg]
%     k1, k2 — OTI 平均值与 5% 分位数的权重
%
%   输出：
%     J — 目标函数值
%         有效圆柱：J = -(k1*mean(OTI) + k2*prctile(OTI,5))
%         无效圆柱：J = 1e6
%
%   说明：
%     本函数使用 20 mm 粗网格进行工作空间搜索，以兼顾 surrogateopt
%     的调用效率。网格设置：
%       seq_x/y = (-400:20:400)*1e-3
%       seq_z   = (-1200:20:-650)*1e-3
%       seq_theta = deg2rad(0:10:30)
%       seq_phi   = deg2rad(-180:30:180)

    params = params(:)';
    H_mm   = params(1);
    h_mm   = params(2);
    r1_mm  = params(3);
    r2_mm  = params(4);
    a_deg  = params(5);
    b_deg  = params(6);

    % ========== 1. 读取并修改参数 ==========
    persistent basic_paras_base
    if isempty(basic_paras_base)
        basic_paras_base = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'm');
    end

    unit_para = basic_paras_base.unit_para;

    % 覆盖待优化参数
    H   = H_mm  * 0.001;
    h   = h_mm  * 0.001;
    r1  = r1_mm * 0.001;
    r2  = r2_mm * 0.001;
    a   = deg2rad(a_deg);
    b   = deg2rad(b_deg);

    l_max = basic_paras_base.l_max;
    l_min = basic_paras_base.l_min;
    R1    = basic_paras_base.R1;
    R2    = basic_paras_base.R2;
    L_tool = basic_paras_base.L_tool;
    joint_u_angle_tilt = basic_paras_base.joint_u_angle_tilt;
    l0_seq = basic_paras_base.l0_seq;

    limb_dir = basic_paras_base.limb_dir;
    limb_dir(4, 2) = deg2rad(90) - a;
    limb_dir(5, 2) = deg2rad(90) + a;
    limb_dir(2, 2) = deg2rad(270) - b;
    limb_dir(3, 2) = deg2rad(270) + b;

    % 重新计算 B, P_m, p_seq
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

    p_seq = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

    % ball screw 方向向量（与 workspace_discrete_v3 一致）
    ball_screw_dir_angle = [deg2rad([35 35 35 35 35]);
                            limb_dir(:, 2)'];
    ball_vector = zeros(3, 5);
    for i = 1 : 5
        ball_vector(1, i) = sin(ball_screw_dir_angle(1,i)) * cos(ball_screw_dir_angle(2,i));
        ball_vector(2, i) = sin(ball_screw_dir_angle(1,i)) * sin(ball_screw_dir_angle(2,i));
        ball_vector(3, i) = cos(ball_screw_dir_angle(1,i));
    end
    ball_vector = -1 * ball_vector;

    % ========== 2. 搜索空间（粗网格） ==========
    seq_x = (-400 : 10 : 400) * unit_para;
    seq_y = (-400 : 10 : 400) * unit_para;
    seq_z = (-1200 : 10 : -650) * unit_para;
    seq_phi   = deg2rad((-180 : 30 : 180));
    seq_theta = deg2rad((0 : 10 : 30));
    len_x = length(seq_x);
    len_y = length(seq_y);
    len_z = length(seq_z);
    len_theta = length(seq_theta);
    len_phi = length(seq_phi);

    ang_threshold = deg2rad(9);

    work_space_ang = [0;0;0;0];
    reachable_thetas = false(len_x, len_y, len_z, len_theta);

    % ========== 3. 工作空间搜索 ==========
    parfor ix = 1 : len_x
        px = seq_x(ix);
        local_ws_ang = [0;0;0;0];
        local_reachable = false(len_y, len_z, len_theta);

        for iy = 1 : len_y
            py = seq_y(iy);
            for iz = 1 : len_z
                pz = seq_z(iz);
                pos_quality = -1;

                for itheta = 1 : len_theta
                    flag_all_phi = 0;
                    theta = seq_theta(itheta);

                    for iphi = 1 : len_phi
                        pos_flag = 0;
                        Pos_ref = [px; py; pz; seq_phi(iphi); theta];
                        T_ref = pos2trans(Pos_ref, B, 'unit', 'rad');
                        vt = T_ref(1:3, 4);
                        R_plant = T_ref(1:3, 1:3);
                        P_v = R_plant * P_m;
                        ball_vector_world = R_plant * ball_vector;

                        for j = 1 : 5
                            vAa = vt + P_v(:, j) - B(:, j);
                            len_vAa = norm(vAa);
                            if (len_vAa >= l_min) && (len_vAa <= l_max)
                                angle_limb_scew = acos(dot(vAa/len_vAa, ball_vector_world(:, j)));
                                if (rad2deg(angle_limb_scew) <= 30) || (j == 1)
                                    pos_flag = pos_flag + 1;
                                end
                            end
                        end

                        if pos_flag < 5
                            flag_all_phi = -1;
                            break;
                        end
                    end  % phi

                    if flag_all_phi < 0
                        break;
                    else
                        pos_quality = theta;
                        local_reachable(iy, iz, itheta) = true;
                    end
                end  % theta

                if pos_quality > ang_threshold
                    local_ws_ang = [local_ws_ang [vt; rad2deg(pos_quality)]];
                end
            end  % z
        end  % y

        work_space_ang = [work_space_ang local_ws_ang(:, 2:end)];
        reachable_thetas(ix, :, :, :) = local_reachable;
    end  % x (parfor)

    % ========== 4. 圆柱提取 ==========
    is_ws_ang = false(len_x, len_y, len_z);
    for k = 2 : size(work_space_ang, 2)
        [~, ix] = min(abs(seq_x - work_space_ang(1, k)));
        [~, iy] = min(abs(seq_y - work_space_ang(2, k)));
        [~, iz] = min(abs(seq_z - work_space_ang(3, k)));
        is_ws_ang(ix, iy, iz) = true;
    end

    best_V = 0; best_R = 0; best_H_cyl = 0;
    best_x0 = 0; best_y0 = 0; best_iz1 = 0; best_iz2 = 0;

    for iz1 = 1 : len_z
        for iz2 = iz1 : len_z
            H_cand = seq_z(iz2) - seq_z(iz1);
            if H_cand < 0.1
                continue;
            end

            common_ws = squeeze(is_ws_ang(:, :, iz1));
            for iz = iz1+1 : iz2
                common_ws = common_ws & squeeze(is_ws_ang(:, :, iz));
            end

            if ~any(common_ws(:))
                continue;
            end

            [reach_ix, reach_iy] = find(common_ws);
            reach_pts = [seq_x(reach_ix).', seq_y(reach_iy).'];

            [unreach_ix, unreach_iy] = find(~common_ws);
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
            if R < 0.055
                continue;
            end
            if (2*R/H_cand) <= 0.6
                continue;
            end

            V = pi * R^2 * H_cand;
            if V > best_V
                best_V = V; best_R = R; best_H_cyl = H_cand;
                best_x0 = best_center_for_interval(1);
                best_y0 = best_center_for_interval(2);
                best_iz1 = iz1; best_iz2 = iz2;
            end
        end
    end

    if best_V <= 0
        J = 1;
        fprintf('[%s] Params=[%.2f,%.2f,%.2f,%.2f,%.2f, %.2f] → No valid cylinder. J=1\n', ...
            string(datetime('now', 'Format', 'HH:mm:ss')), H_mm, h_mm, r1_mm, r2_mm, a_deg, b_deg);
        return;
    end

    R_cyl = best_R;
    z_min_cyl = seq_z(best_iz1);
    z_max_cyl = seq_z(best_iz2);

    % ========== 5. OTI 计算 ==========
    % 收集圆柱内所有需要评估的 [x,y,z,theta,phi]
    oti_all = [];

    parfor ix = 1 : len_x
        px = seq_x(ix);
        local_otis = [];
        for iy = 1 : len_y
            py = seq_y(iy);
            if sqrt((px-best_x0)^2 + (py-best_y0)^2) > R_cyl
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

                    oti_phi = zeros(len_phi, 1);
                    for iphi = 1 : len_phi
                        Pos_ref = [px; py; pz; seq_phi(iphi); theta];
                        T_ref = pos2trans(Pos_ref, B, 'unit', 'rad');
                        [~, ~, oti_val, ~] = compute_ltigci(T_ref, B, l0_seq, P_m, p_seq);
                        oti_phi(iphi) = oti_val;
                    end
                    local_otis = [local_otis; mean(oti_phi)];
                end
            end
        end
        oti_all = [oti_all; local_otis];
    end

    if isempty(oti_all)
        J = 1e6;
        fprintf('[%s] Params=[%.2f,%.2f,%.2f,%.2f,%.2f,%.2f] → Cylinder found but OTI empty. J=1e6\n', ...
            string(datetime('now', 'Format', 'HH:mm:ss')), H_mm, h_mm, r1_mm, r2_mm, a_deg, b_deg);
        return;
    end

    mean_otis = mean(oti_all);
    prc5_otis = prctile(oti_all, 5);
    J = -(k1 * mean_otis + k2 * prc5_otis);

    fprintf('[%s] Params=[%.2f,%.2f,%.2f,%.2f,%.2f,%.2f] | Cyl:R=%.4f,H=%.4f | OTI_mean=%.4f,OTI_5%%=%.4f | J=%.6f\n', ...
        string(datetime('now', 'Format', 'HH:mm:ss')), H_mm, h_mm, r1_mm, r2_mm, a_deg, b_deg, ...
        R_cyl, best_H_cyl, mean_otis, prc5_otis, J);

end
