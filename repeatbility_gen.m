% 生成重复性测试用轨迹
% 共有5个点，以相同的方式一次经过，循环30次

clear
path_add()

flag = 1;  % 1表示平面空间，2表示球面空间


if flag == 1
    a_square = 80;
    Pc = [0; 0; -972.5];
else
    a_square = 45;
    Pc = [0; -190; -877.5];
end


P1 = [0; 0; 0];
P2 = [0.4*a_square; 0.4*sqrt(2)*a_square; 0.4*sqrt(2)*a_square];
P3 = [-0.4*a_square; 0.4*sqrt(2)*a_square; 0.4*sqrt(2)*a_square];
P4 = [-0.4*a_square; -0.4*sqrt(2)*a_square; -0.4*sqrt(2)*a_square];
P5 = [0.4*a_square; -0.4*sqrt(2)*a_square; -0.4*sqrt(2)*a_square];
P = [P1 P2 P3 P4 P5] + Pc;

%% 1. 添加姿态信息并循环30次，生成刀具末端点位姿序列 (x, y, z, phi, theta)
% 每个点的姿态均为 phi = 0, theta = 0；角度单位：deg
n_repeat = 30;
tool_points = [];  % 每列 [x; y; z; phi; theta]
for k = 1 : n_repeat
    for i = 1 : size(P, 2)
        tool_points = [tool_points [P(:, i); 0; 0]];
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
    write_cmd_csv(fullfile('calibration_data', 'repeat_tool_points_planar.csv'), tool_points);
    write_cmd_csv(fullfile('calibration_data', 'repeat_surf_points_planar.csv'), surf_points);
elseif flag == 2
    write_cmd_csv(fullfile('calibration_data', 'repeat_tool_points_spherical.csv'), tool_points);
    write_cmd_csv(fullfile('calibration_data', 'repeat_surf_points_spherical.csv'), surf_points);
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
    % pts: 每列 [x; y; z; phi; theta]，长度单位mm；输出csv时转换为um
    fid = fopen(file_name, 'w');
    fprintf(fid, 'cmd,x,y,z,phi,theta,ticks\n');
    for i = 1 : size(pts, 2)
        fprintf(fid, '2,%g,%g,%g,%g,%g,0\n', pts(1:3, i)*1000, pts(4, i), pts(5, i));
        fprintf(fid, '20,0,0,0,0,0,5000\n');  % 中间暂停
    end
    fprintf(fid, '2,0,0,-800000,0,0,0\n');
    fprintf(fid, '0,0,0,0,0,0,0\n');  % 用来表示结束
    fclose(fid);
end





