clear
fig_rotation_show = 0;  % 1开启展示旋转
gif_generate_flag = 0;  % 1为开启录制功能，运行一次程序后记得改文件名
flag_range_plot = 0;  % 1为开启，绘制出指定高度的空间范围

path_add()
fprintf('>>>= start (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));
%--------parameter3--------
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

L_tool = 0;

pos_plant = [0; 0; -800]*unit_para;  % 后面作图用，不参与空间搜索
alpha_plant = 0;  % 绕 x
beta_plant = 0;  % 绕 y
gamma_plant = 0;  % 绕 z


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

P_v = zeros(3, 5);  % 只变换了方向，没变换起点
P = zeros(3, 5);    % 末端点坐标
for i = 1 : 5
    P_v(:, i) = R_plant * P_m(:, i);
    P(:, i) = P_v(:, i) + pos_plant;
end


% ball screw dir vector
% ball_screw_dir_angle_deg = [145 145 145  145 145;   % 与Z轴夹角  % pi - 35°
%                             -90  30 150 -150 -30];  % 与x轴夹角  % (dir+180)mod360
% ball_screw_dir_angle = ball_screw_dir_angle_deg / 180 * pi;

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

% -----end-parameter3------



% ------search space-------
seq_x = (-400 : 10 : 400)*unit_para;
seq_y = (-400 : 10 : 400)*unit_para;
seq_z = (-1200 : 10 : -650)*unit_para;
seq_phi = deg2rad((-180 : 30 : 180));
seq_theta = deg2rad((0 : 2 : 30));
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


%% 遍历
parfor ix = 1 : len_x
    px = seq_x(ix);
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
                    R_plant = T_ref(1:3, 1:3);  % 由于考虑了支链约束，因此工作空间计算与上一版不同
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


                    % 暂时没用到
                    % s_limb_move = zeros(3, 5);  % 动平台坐标系下支链方向向量
                    % for i_limb = 1 : 5
                    %     s_limb_move(:, i_limb) = R_plant' * s_limb(:, i_limb);
                    % end


                    % % 雅克比矩阵条件数计算
                    % x_m = [1; 0; 0];
                    % J1 = [s_limb (R_plant * x_m)];
                    % J2 = [cross(P_v, s_limb, 1)  ...
                    % (cross(P_v(:, 1 ), (R_plant * x_m)) + l_limb(1)*cross((R_plant * x_m), s_limb(:, 1)))];
                    % J = [J1' J2'];
                    % % structure of J
                    % % J1'  J2'
                    % % s1   \arr {ObP1} \times s1
                    % % s2   \arr {ObP2} \times s2
                    % % s3   \arr {ObP3} \times s3
                    % % s4   \arr {ObP4} \times s4
                    % % s5   \arr {ObP5} \times s5
                    % % xb   L1(xb \times s1) + \arr {ObP1} \times xb
                    % cond_J = cond(J);
                    % det_J = det(J);
                    % if abs(det_J) < 1e-8
                    %     fprintf("x=%d,y=%d,z=%d,det_J=%.4f\n",seq_x(ix),seq_y(iy),seq_z(iz),det_J);
                    % end
            
                end  % phi

                % 记录theta情况
                if(flag_all_phi < 0)
                    break;
                else
                    pos_quality = theta;
                end
            end  % theta

            % 含角度的工作空间
            if pos_quality > ang_threshold
                p_mark = [vt; rad2deg(pos_quality)];
                work_space_ang = [work_space_ang p_mark];
            end


            % 工作空间轮廓：只画z方向的上下两端点
            if pos_quality > -1  % 如果所有条件均允许
                if z_min_point == zeros(3,1)  % 如果第一次进入循环
                    z_min_point = vt;
                    z_max_point = vt;
                else
                    z_max_point = vt;
                end

                pos_count = pos_count + 1;
                % work_space = [ work_space_ang vt];  % 全局搜索
            end


        end  % z


        % 工作空间轮廓：将工作空间界限添加到最后的作图中
        if (z_min_point(3) ~= 0) && (z_max_point(3) ~= 0)
            work_space_up = [work_space_up z_max_point];
            work_space_down = [work_space_down z_min_point];
        end




    end  % y
end  % x

%% ------- plot -------
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

%% 含角度的工作空间
fig_ang = figure('Color', [1 1 1]);
scatter3(work_space_ang(1,2:end), work_space_ang(2,2:end), work_space_ang(3,2:end), 2, work_space_ang(4,2:end), 'filled');
grid on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
colormap(parula); % 可选：jet, parula, turbo, hot等
colorbar;

%% range plot
if flag_range_plot == 1
    fig2 = figure('Color', [1 1 1]);
    n = length(work_space_down(1,:));
    for i = 1 : n
        if work_space_down(3, i) == -820*unit_para
            plot(work_space_down(1,i), work_space_down(2,i),'*')
            hold on
        end
    end
    for i = 1 : n
        if work_space_down(3, i) == -850*unit_para
            plot(work_space_down(1,i), work_space_down(2,i),'*')
            hold on
        end
    end
    axis equal
end

%% view rotation
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


fprintf('>>>= done (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));



