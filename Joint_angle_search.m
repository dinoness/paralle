clear
fprintf('>>>= Joint_angle_search Start (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));
workspace_plot_flag = 0;  % 1开启工作区间绘制
fig_rotation_show = 0;  % 1开启展示旋转
gif_generate_flag = 0;  % 1为开启录制功能，运行一次程序后记得改文件名


%--------parameter3--------
T = readtable('parameters.xlsx', 'Range', 'A2:B12');
paras = table2array(T(:, 2));
l_max = paras(1);
l_min = paras(2);  % 670
R1 = paras(3);  % 800
R2 = paras(4);  % 600
H = paras(5);  % 20
r1 = paras(6);  % 100
r2 = paras(7);  % 80
h = paras(8);  % 100

pos_plant = [0; 0; -600];  % 后面作图用，不参与空间搜索
alpha_plant_deg = paras(9);  % 绕 x
beta_plant_deg = paras(10);  % 绕 y
gamma_plant_deg = paras(11);  % 绕 z
alpha_plant = alpha_plant_deg / 180 * pi;  % 绕 x
beta_plant = beta_plant_deg / 180 * pi;  % 绕 y
gamma_plant = gamma_plant_deg / 180 * pi;  % 绕 z

% l_max = 1030;
% l_min = 720;  % 670
% R1 = 700;  % 800
% R2 = 600;  % 600
% H = -60;  % 20
% r1 = 120;  % 100
% r2 = 100;  % 80
% h = 20;  % 100

% static plant
B1 = [R1*cos(pi/2);   R1*sin(pi/2);   0];
B2 = [R1*cos(7*pi/6); R1*sin(7*pi/6); 0];
B3 = [R1*cos(-pi/6);  R1*sin(-pi/6);  0];
B4 = [R2*cos(pi/6);   R2*sin(pi/6);   H];
B5 = [R2*cos(5*pi/6); R2*sin(5*pi/6); H];
B = [B1 B2 B3 B4 B5];


% Position and Posture of Move Plant
% pos_plant = [0; 0; -600];  % 后面作图用，不参与空间搜索
% alpha_plant = 0 / 180 * pi;  % 绕 x
% beta_plant = 0 / 180 * pi;  % 绕 y
% gamma_plant = 0 / 180 * pi;  % 绕 z

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

% move plant parameter
P1_m = [r1*cos(pi/2);   r1*sin(pi/2);   0];
P2_m = [r1*cos(7*pi/6); r1*sin(7*pi/6); 0];
P3_m = [r1*cos(-pi/6);  r1*sin(-pi/6);  0];
P4_m = [r2*cos(pi/6);   r2*sin(pi/6);   h];
P5_m = [r2*cos(5*pi/6); r2*sin(5*pi/6); h];
P_m = [P1_m P2_m P3_m P4_m P5_m];
P_v = zeros(3, 5);  % 只变换了方向，没变换起点
P = zeros(3, 5);    % 末端点坐标
for i = 1 : 5
    P_v(:, i) = R_plant * P_m(:, i);
    P(:, i) = P_v(:, i) + pos_plant;
end

% ball screw dir vector
ball_screw_dir_angle_deg = [145 145 145  145 145; 
                            -90  30 150 -150 -30];
ball_screw_dir_angle = ball_screw_dir_angle_deg / 180 * pi;
ball_vector = zeros(3, 5);
ball_vector_world = zeros(3, 5);

static_joint_dir_angle_deg = [ 0  0  0  0  0;
                              -90  30 150 -150 -30];
static_joint_dir_angle = static_joint_dir_angle_deg / 180 * pi;
static_joint_vector = zeros(3, 5);

for i_ball = 1 : 5
    ball_vector(1, i_ball) = sin(ball_screw_dir_angle(1, i_ball)) * cos(ball_screw_dir_angle(2, i_ball));
    ball_vector(2, i_ball) = sin(ball_screw_dir_angle(1, i_ball)) * sin(ball_screw_dir_angle(2, i_ball));
    ball_vector(3, i_ball) = cos(ball_screw_dir_angle(1, i_ball));
    ball_vector_world = R_plant * ball_vector;

    static_joint_vector(1, i_ball) = sin(static_joint_dir_angle(1, i_ball)) * cos(static_joint_dir_angle(2, i_ball));
    static_joint_vector(2, i_ball) = sin(static_joint_dir_angle(1, i_ball)) * sin(static_joint_dir_angle(2, i_ball));
    static_joint_vector(3, i_ball) = cos(static_joint_dir_angle(1, i_ball));
end

% -----end-parameter3------



% ------search space-------
% 搜索的是核心工作区间，不一定是完整的工作区间
seq_x = -200 : 10 : 200;
seq_y = -200 : 10 : 200;
seq_z = -900 : 10 : -700;

% assistant parameter
wors_space = [];
work_space_up = [0;0;0];
work_space_down = [0;0;0];
limb_angle = [];
limb_angle_2_static = [];
s_limb_seq = cell(1, 5);
for i = 1 : 5
    s_limb_seq{i} = [];
end


% length transform
for ix = 1 : length(seq_x)
    for iy = 1 : length(seq_y)
        % x,y方向进行遍历

        z_min_point = zeros(3,1);
        z_max_point = zeros(3,1);

        % 搜索z向上的工作空间界限
        for iz = 1 : length(seq_z)
            vt = [seq_x(ix); seq_y(iy); seq_z(iz)];  % 搜索的目标点
            pos_flag = 0;  % 位置可达标志位
            s_limb = zeros(3, 5);  % 支链的方向向量
            
            for j = 1 : length(P_v(1, :))
                vAa = vt + P_v(:, j) - B(:, j);
                len_vAa = norm(vAa);
                s_limb(:, j) = vAa / len_vAa;


                if (len_vAa >= l_min)&&(len_vAa <= l_max)
                    pos_flag = pos_flag + 1;
                end
            end


            % 只画z方向的上下两端点
            if pos_flag == length(P_v(1, :))  % 如果所有条件均允许
                
                % 只计算在工作空间中的点
                s_limb_move = zeros(3, 5);
                limb_angle_in = zeros(1, 5);
                limb_angle_2_static_in = zeros(1, 5);
                for i_limb = 1 : 5
                    s_limb_move_in = R_plant' * s_limb(:, i_limb);
                    limb_angle_in(1, i_limb) = acos(dot(s_limb_move_in, ball_vector_world(:, i_limb)));
                    limb_angle_2_static_in(1, i_limb) = acos(dot(s_limb_move_in, static_joint_vector(:, i_limb)));
                    s_limb_move(:, i_limb) = s_limb_move_in;
                    s_limb_seq{i_limb} = [s_limb_seq{i_limb} s_limb_move(:, i_limb)];
                end

                limb_angle = [limb_angle ; limb_angle_in];
                limb_angle_2_static = [limb_angle_2_static; limb_angle_2_static_in];


                if z_min_point == zeros(3,1)  % 如果第一次进入循环
                    z_min_point = vt;
                    z_max_point = vt;
                else
                    z_max_point = vt;
                end
                % work_space = [wors_space vt];  % 全局搜索
            end
      
        end

        % 将工作空间界限添加到最后的作图中
        if workspace_plot_flag == 1
            if (z_min_point(3) ~= 0) && (z_max_point(3) ~= 0)
                work_space_up = [work_space_up z_max_point];
                work_space_down = [work_space_down z_min_point];
            end
        end
    end
end


% -------- angle data process -----------
limb_angle_deg = limb_angle / pi * 180;
limb_angle_mean = mean(limb_angle_deg);
limb_angle_var = var(limb_angle_deg);
limb_angle_max = max(limb_angle_deg);
limb_angle_min = min(limb_angle_deg);
limb_2_static_angle_deg = limb_angle_2_static / pi * 180;
limb_2_static_angle_max = max(limb_2_static_angle_deg);
limb_2_static_angle_min = min(limb_2_static_angle_deg);

num_point = length(s_limb_seq{1}(1, :));
value_vec = zeros(5, num_point);  % 每个向量与其他向量点积之和
idx_vec = zeros(1, 5);
desired_vec = zeros(3, 5);
desired_ball_screw_dir_angle = zeros(2, 5);  % phi; theta
for i_limb = 1 : 5
    s_limb_i_seq = s_limb_seq{i_limb};

    for i_vector = 1 : num_point
        v_limb = s_limb_i_seq(:, i_vector);
        for i_vector2 = 1 : num_point
            value_vec(i_limb, i_vector) = value_vec(i_limb, i_vector) + dot(v_limb, s_limb_i_seq(:, i_vector2));
        end
    end

    [max_val, idx_vec(i_limb)] = max(value_vec(i_limb, :));
    desired_vec(:, i_limb) = s_limb_i_seq(:, idx_vec(i_limb));

    % 用球面坐标系来描述方位向量
    phi = acos(desired_vec(3, i_limb));
    theta = acos(desired_vec(1, i_limb) / sin(phi));
    desired_ball_screw_dir_angle(1, i_limb) = phi;
    desired_ball_screw_dir_angle(2, i_limb) = theta;
end

desired_ball_screw_dir_angle_deg = desired_ball_screw_dir_angle / pi * 180;
fprintf("rotation of plate: x - %.2f°, y - %.2f°, z - %.2f°\n",alpha_plant_deg ,beta_plant_deg, gamma_plant_deg);
fprintf("phi(angle to z axis), theta(angle to x asix in xoy plant)\n");
for i_limb = 1 : 5
    fprintf("limb %d, ", i_limb);
    fprintf("phi = %.4f, ", desired_ball_screw_dir_angle_deg(1, i_limb));
    fprintf("theta = %.4f, ", desired_ball_screw_dir_angle_deg(2, i_limb));
    fprintf("max_angle = %.4f, ", limb_angle_max(i_limb));
    fprintf("min_angle = %.4f\n", limb_angle_min(i_limb));
end


% ---------------- plot -----------------
if workspace_plot_flag == 1
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
    % 工作空间散点
    % plot3(work_space_up(1,2:end), work_space_up(2,2:end), work_space_up(3,2:end), '.', 'Color', '#00BFFF');
    % plot3(work_space_down(1,2:end), work_space_down(2,2:end), work_space_down(3,2:end), '.', 'Color', '#4169E1');
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
