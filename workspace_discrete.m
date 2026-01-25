clear
fig_rotation_show = 0;  % 1开启展示旋转
gif_generate_flag = 0;  % 1为开启录制功能，运行一次程序后记得改文件名

% 还没有考虑关节角度的限制
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

pos_plant = [0; 0; -800];  % 后面作图用，不参与空间搜索
alpha_plant = paras(9) / 180 * pi;  % 绕 x
beta_plant = paras(10) / 180 * pi;  % 绕 y
gamma_plant = paras(11) / 180 * pi;  % 绕 z

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
ball_screw_dir_angle_deg = [145 145 145  145 145;   % 与Z轴夹角
                            -90  30 150 -150 -30];  % 与x轴夹角
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
seq_x = -400 : 10 : 400;
seq_y = -400 : 10 : 400;
seq_z = -1000 : 10 : -500;

% assistant parameter
wors_space = [];
work_space_up = [0;0;0];
work_space_down = [0;0;0];
pos_count = 0;  % 空间点计数

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
            l_limb = zeros(1, 5);  % 支链长度
            
            for j = 1 : length(P_v(1, :))
                vAa = vt + P_v(:, j) - B(:, j);
                len_vAa = norm(vAa);
                l_limb(j) = len_vAa;
                s_limb(:, j) = vAa / len_vAa;

                if (len_vAa >= l_min)&&(len_vAa <= l_max)  % ===========支链长度条件===========
                    angle_limb_scew = acos(dot(s_limb(:, j), ball_vector_world(:, j)));  % 支链与球铰轴线夹角
                    angle_limb_scew_deg = angle_limb_scew / pi * 180;
                    if(angle_limb_scew_deg <= 30 )  % ===========关节角度条件===========
                        pos_flag = pos_flag + 1;
                    end
                end
            end

            s_limb_move = zeros(3, 5);  % 动平台坐标系下支链方向向量
            for i_limb = 1 : 5
                s_limb_move(:, i_limb) = R_plant' * s_limb(:, i_limb);
            end


            % 只画z方向的上下两端点
            
            if pos_flag == length(P_v(1, :))  % 如果所有条件均允许
                if z_min_point == zeros(3,1)  % 如果第一次进入循环
                    z_min_point = vt;
                    z_max_point = vt;
                else
                    z_max_point = vt;
                end

                pos_count = pos_count + 1;
                % work_space = [wors_space vt];  % 全局搜索
            end

            % 雅克比矩阵条件数计算
            x_m = [1; 0; 0];
            J1 = [s_limb (R_plant * x_m)];
            J2 = [cross(P_v, s_limb, 1)  ...
            (cross(P_v(:, 1 ), (R_plant * x_m)) + l_limb(1)*cross((R_plant * x_m), s_limb(:, 1)))];
            J = [J1' J2'];
            % structure of J
            % J1'  J2'
            % s1   \arr {ObP1} \times s1
            % s2   \arr {ObP2} \times s2
            % s3   \arr {ObP3} \times s3
            % s4   \arr {ObP4} \times s4
            % s5   \arr {ObP5} \times s5
            % xb   L1(xb \times s1) + \arr {ObP1} \times xb
            cond_J = cond(J);
            det_J = det(J);
            if abs(det_J) < 0.01
                fprintf("x=%d,y=%d,z=%d,det_J=%.4f\n",seq_x(ix),seq_y(iy),seq_z(iz),det_J);
            end
      
        end

        % 将工作空间界限添加到最后的作图中
        if (z_min_point(3) ~= 0) && (z_max_point(3) ~= 0)
            work_space_up = [work_space_up z_max_point];
            work_space_down = [work_space_down z_min_point];
        end
    end
end

% ------- plot -------
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


fprintf('>>>= workspace_discrete done (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));



