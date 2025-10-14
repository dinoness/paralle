clear
workspace_plot_flag = 0;  % 1开启工作区间绘制
fig_rotation_show = 0;  % 1开启展示旋转
gif_generate_flag = 0;  % 1为开启录制功能，运行一次程序后记得改文件名


%--------parameter3--------
l_max = 1000;
l_min = 670;
R1 = 800;  % 800
R2 = 600;  % 600
H = 20;  % 20
r1 = 100;
r2 = 80;  % 80
h = 100;  % 100

% static plant
B1 = [R1*cos(pi/2);   R1*sin(pi/2);   0];
B2 = [R1*cos(7*pi/6); R1*sin(7*pi/6); 0];
B3 = [R1*cos(-pi/6);  R1*sin(-pi/6);  0];
B4 = [R2*cos(pi/6);   R2*sin(pi/6);   H];
B5 = [R2*cos(5*pi/6); R2*sin(5*pi/6); H];
B = [B1 B2 B3 B4 B5];


% Position and Posture of Move Plant
pos_plant = [0; 0; -600];  % 后面作图用，不参与空间搜索
alpha_plant = 0 / 180 * pi;  % 绕 x
beta_plant = 0 / 180 * pi;  % 绕 y
gamma_plant = 0 / 180 * pi;  % 绕 z

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
% -----end-parameter3------



% ------search space-------
seq_x = -200 : 10 : 200;
seq_y = -200 : 10 : 200;
seq_z = -650 : 10 : -400;

% assistant parameter
wors_space = [];
work_space_up = [0;0;0];
work_space_down = [0;0;0];
limb_angle = [];
s_limb_2_seq = [];


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
                for i_limb = 1 : 5
                    s_limb_move_in = R_plant' * s_limb(:, i_limb);
                    limb_angle_in(1, i_limb) = acos(dot(s_limb_move_in, [0; 0; -1]));
                    s_limb_move(:, i_limb) = s_limb_move_in;
                end

                limb_angle = [limb_angle ; limb_angle_in];
                s_limb_2_seq = [s_limb_2_seq s_limb_move(:, 2)];


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

num_point = length(s_limb_2_seq(1, :));
value_vec = zeros(1, num_point);
for i_vector = 1 : num_point
    v_limb = s_limb_2_seq(:, i_vector);
    for i_vector2 = 1 : num_point
        value_vec(i_vector) = value_vec(i_vector) + dot(v_limb, s_limb_2_seq(:, i_vector2));
    end
end
[max_val, idx] = max(value_vec)
s_limb_2_seq(:, idx)

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

fprintf('>>>= Joint_angle_search done =<<<\n');



