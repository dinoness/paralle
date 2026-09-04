% 重复性计算脚本（ISO 9283 位姿重复性 / 位置重复性 RP）
% 相同路径经过 5 个点，重复 30 次；数据文件 repeatability.txt
% 首行与末行为参考点（同一物理点的复测），去除后剩余 150 行 = 30 循环 × 5 点
%
% ISO 9283 位置重复性定义：
%   l_i = |p_i - p̄|          （各次测量到集群中心的距离）
%   RP  = l̄ + 3 * S_l        （l̄ 为 l_i 均值，S_l 为其标准差）
clear

%% 读取数据
file_name = fullfile('repeatability', 'repeatability.txt');
fid = fopen(file_name, 'r');
if fid < 0
    error('无法打开文件: %s', file_name);
end
pts = [];
while true
    ln = fgetl(fid);
    if ~ischar(ln)
        break;
    end
    if isempty(strtrim(ln))
        continue;
    end
    tokens = regexp(ln, ',', 'split');
    pts(:, end+1) = str2double(strtrim(tokens(2:4))).'; %#ok<AGROW>
end
fclose(fid);

% 去除首行与末行（参考点复测）
pts = pts(:, 2:end-1);

n_point = 5;
n_repeat = 30;
if size(pts, 2) ~= n_point * n_repeat
    error('数据点数 %d 与预期 %d 不符', size(pts, 2), n_point * n_repeat);
end

%% 按点分组：每循环依次经过 P1~P5
% groups{j}: 3×n_repeat，第 j 个点的 30 次测量
groups = cell(1, n_point);
for j = 1 : n_point
    groups{j} = pts(:, j : n_point : end);
end

%% ISO 9283 位置重复性计算
fprintf('%-6s %12s %12s %12s %12s %12s %12s\n', ...
    '点', 'l_bar(mm)', 'S_l(mm)', 'RP(mm)', 'S_x(mm)', 'S_y(mm)', 'S_z(mm)');
RP_seq = zeros(1, n_point);
for j = 1 : n_point
    p = groups{j};
    p_bar = mean(p, 2);                    % 集群中心
    l = sqrt(sum((p - p_bar).^2, 1));      % 各次测量到中心的距离
    l_bar = mean(l);
    S_l = std(l);
    RP = l_bar + 3 * S_l;
    RP_seq(j) = RP;
    S_xyz = std(p, 0, 2);                  % 各坐标分量标准差
    fprintf('P%-5d %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
        j, l_bar, S_l, RP, S_xyz(1), S_xyz(2), S_xyz(3));
end
fprintf('最大位置重复性 RP_max = %.4f mm\n', max(RP_seq));

%% 可视化：各点 30 次测量的散布
figure('Color', [1 1 1]);
hold on
colors = lines(n_point);
for j = 1 : n_point
    p = groups{j};
    p_bar = mean(p, 2);
    dp = (p - p_bar) * 1000;  % 相对集群中心的偏差，um
    plot3(dp(1,:), dp(2,:), dp(3,:), '.', 'Color', colors(j,:), 'MarkerSize', 10);
end
grid on
axis equal
xlabel('\Deltax (um)')
ylabel('\Deltay (um)')
zlabel('\Deltaz (um)')
title('各测试点相对集群中心的散布')
legend({'P1','P2','P3','P4','P5'}, 'Location', 'best')

fprintf('>>>= done (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));
