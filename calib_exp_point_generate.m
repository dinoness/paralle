clear
path_add()

%% 生成刀具空间
flag = 2;  % 1表示平面空间，2表示球面空间

% step_size = 15;
% x_seq = -105 : step_size : 105;
% y_seq = -105 : step_size : 105;

if flag == 1
    D = 200;
    R = D / 2;
    H = 15;
    step_size = 25;
    x_seq = -200 : step_size : 200;
    y_seq = -200 : step_size : 200;
    z_seq = -40 : 20 : 40;
    x0 = 0;
    y0 = 0;
    z0 = -972.5;
elseif flag == 2
    D = 100;
    R = D / 2;
    H = 15;
    step_size = 20;
    x_seq = -200 : step_size : 200;
    y_seq = -200 : step_size : 200;
    z_seq = -40 : step_size : 40;
    x0 = 0;
    y0 = 0;
    z0 = -877.5;
end

x_length = length(x_seq);
y_length = length(y_seq);
z_length = length(z_seq);
points = [];

for iz = 1 : z_length
    for ix = 1 : x_length
        for iy = 1 : y_length
            p = [x_seq(ix)+x0; y_seq(iy)+y0; z_seq(iz)+z0];
            dis = norm(p - [x0; y0; z0]);  % 到空间中心 (x0,y0,z0) 的距离
            if dis <= R
                points = [points p];
            end
        end
    end
end
plot(points(1,:),points(2,:),'.');
fprintf('point num = %d \n', length(points(1,:)))

%% 1. 添加姿态信息，生成刀具末端点位姿序列 (x, y, z, phi, theta)
% 平面空间：phi = 0, theta = 0
% 球面空间：theta = 0 取一个姿态；theta = 5°, 10° 时 phi 分别取 -120°, 0°, 120°
% 角度单位：deg
if flag == 1
    pose_seq = [0 0];  % [phi theta]
elseif flag == 2
    pose_seq = [   0  0;
                -120  5;
                   0  5;
                 120  5;
                -120 10;
                   0 10;
                 120 10];
end

tool_points = [];  % 每列 [x; y; z; phi; theta]
for i = 1 : size(points, 2)
    for j = 1 : size(pose_seq, 1)
        tool_points = [tool_points [points(:, i); pose_seq(j, 1); pose_seq(j, 2)]];
    end
end
fprintf('tool point num (with pose) = %d \n', size(tool_points, 2))

%% 2. 刀具末端点序列 -> 动平台表面中心点序列
% 动平台坐标系约定（单位mm）：原点在动平台中心，
% 刀具末端点坐标 (0,0,-150)，动平台表面中心点坐标 (0,0,-7.5)
% 位姿转移矩阵参考 pos2trans：z轴由 phi/theta 确定，x轴由 B1 与原点连线方向确定
basic_paras = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'mm');
B = basic_paras.B;

L_tip2center = 150;   % 刀尖到动平台中心的距离
L_center2surf = 7.5;  % 动平台中心到表面中心的距离

surf_points = [];  % 每列 [x; y; z; phi; theta]，phi/theta 与刀尖位姿一致
for i = 1 : size(tool_points, 2)
    p_tip = tool_points(1:3, i);
    phi = tool_points(4, i);
    theta = tool_points(5, i);

    % 动平台中心 = 刀尖沿刀具轴线（z轴）上移 150 mm
    z_axis = [sind(theta)*cosd(phi); sind(theta)*sind(phi); cosd(theta)];
    t_center = p_tip + L_tip2center * z_axis;

    % 以动平台中心为原点构建位姿转移矩阵，求表面中心点的世界坐标
    T = pos2trans([t_center; phi; theta], B, 'unit', 'deg');
    p_surf = T(1:3, 1:3) * [0; 0; -L_center2surf] + T(1:3, 4);

    surf_points = [surf_points [p_surf; phi; theta]];
end
fprintf('surface center point num = %d \n', size(surf_points, 2))

%% 输出与可视化
% csv格式：cmd,x,y,z,phi,theta,ticks
% 每个点：cmd=2, ticks=0；点与点之间插入一行 20,0,0,0,0,0,5000；
% 最后一行为 2,0,0,-800,0,0,0
if ~exist('calibration_data', 'dir')
    mkdir('calibration_data');
end
if flag == 1
    write_cmd_csv(fullfile('calibration_data', 'calib_exp_tool_points_planar.csv'), tool_points);
    write_cmd_csv(fullfile('calibration_data', 'calib_exp_surf_points_planar.csv'), surf_points);
elseif flag == 2
    write_cmd_csv(fullfile('calibration_data', 'calib_exp_tool_points_spherical.csv'), tool_points);
    write_cmd_csv(fullfile('calibration_data', 'calib_exp_surf_points_spherical.csv'), surf_points);
end

figure('Color', [1 1 1]);
plot3(surf_points(1,:), surf_points(2,:), surf_points(3,:), '-');
grid on
axis equal
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
title('Surface center points')

fprintf('>>>= done (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));

function write_cmd_csv(file_name, pts)
    % pts: 每列 [x; y; z; phi; theta]
    fid = fopen(file_name, 'w');
    fprintf(fid, 'cmd,x,y,z,phi,theta,ticks\n');
    for i = 1 : size(pts, 2)
        fprintf(fid, '2,%g,%g,%g,%g,%g,0\n', pts(:, i));
        if i < size(pts, 2)
            fprintf(fid, '20,0,0,0,0,0,5000\n');
        end
    end
    fprintf(fid, '2,0,0,-800,0,0,0\n');
    fclose(fid);
end
