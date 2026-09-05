function [pose_seq, T_seq] = calib_pts2pose_seq(data_dir, base_file)
%CALIB_PTS2POSE_SEQ 处理标定测量数据，生成动平台表面坐标系位姿序列
%   输入：
%     data_dir  — 数据目录，包含 t1.txt / t2.txt / t3.txt 三个文件，
%                 每行格式：点编号, x, y, z（默认 'calib_p'）
%     base_file — 世界坐标系定义文件（默认 'base.txt'），包含三个点：
%                 三点平均值为原点，p2指向p3的方向为x轴，
%                 p1到x轴的垂线（由垂足指向p1）为y轴，z轴由右手定则确定
%   输出：
%     pose_seq — 5×n 位姿序列，每列 [x; y; z; phi; theta]，
%                为动平台表面坐标系在世界坐标系下的位姿，
%                长度单位 mm，角度单位 deg
%     T_seq    — 4×4×n 齐次变换矩阵序列，由测量点构建的坐标系直接
%                生成，为动平台表面坐标系在世界坐标系下的表示，
%                平移单位 mm
%
%   处理流程：
%     1. 读取 t1/t2/t3 测量文件，去除首行与末行（参考点复测），
%        其余行为测量点；
%     2. 由 base.txt 的三个点构建世界坐标系；
%     3. 相同编号的三个测量点构建测量子坐标系：原点为三点均值，
%        x轴为 p3 指向 p2 的方向，y轴为 p1 到x轴的垂线方向（由p1指向
%        垂足），z轴由右手定则确定；
%     4. 动平台表面坐标系：坐标轴与测量子坐标系平行，其原点在测量子坐标
%        系中的位置为 (0, 0, 37.8) mm；最终输出动平台表面坐标系在世界
%        坐标系下的位姿序列（位置和姿态矩阵均在世界系下描述）。

    if nargin < 1 || isempty(data_dir)
        data_dir = 'calib_p';
    end
    if nargin < 2 || isempty(base_file)
        base_file = 'base.txt';
    end

    O_surf_in_meas = [0; 0; 37.8];  % 动平台表面原点在测量子坐标系中的位置 (mm)

    % 1. 测量点（去除首末行）
    [ids1, pts1] = read_points_file(fullfile(data_dir, 't1.txt'));
    [ids2, pts2] = read_points_file(fullfile(data_dir, 't2.txt'));
    [ids3, pts3] = read_points_file(fullfile(data_dir, 't3.txt'));

    ids1 = ids1(2:end-1);  pts1 = pts1(:, 2:end-1);
    ids2 = ids2(2:end-1);  pts2 = pts2(:, 2:end-1);
    ids3 = ids3(2:end-1);  pts3 = pts3(:, 2:end-1);

    % 2. 世界坐标系（base_file 三点定义）
    [~, base_pts] = read_points_file(fullfile(data_dir, base_file));
    if size(base_pts, 2) ~= 3
        error('%s 应包含 3 个点，实际读取到 %d 个', base_file, size(base_pts, 2));
    end
    [R_w, O_w] = build_world_frame(base_pts(:,1), base_pts(:,2), base_pts(:,3));

    % 3. 相同编号的测量点构建测量子坐标系，转换到世界坐标系下
    common_ids = intersect(intersect(ids1, ids2), ids3);
    n = numel(common_ids);
    pose_seq = zeros(5, n);
    T_seq = zeros(4, 4, n);
    for k = 1 : n
        id = common_ids(k);
        q1 = pts1(:, ids1 == id);
        q2 = pts2(:, ids2 == id);
        q3 = pts3(:, ids3 == id);

        [R_sub, O_sub] = build_frame(q1, q2, q3);

        % 测量子坐标系在世界坐标系下的表示
        R_rel = R_w.' * R_sub;
        t_rel = R_w.' * (O_sub - O_w);

        % 4. 测量子坐标系 -> 动平台表面坐标系（纯平移，坐标轴平行）
        t_world = t_rel + R_rel * O_surf_in_meas;
        T_seq(:, :, k) = [R_rel, t_world; 0 0 0 1];

        % 由旋转矩阵的 z 轴列恢复 phi/theta（与 pos2trans 约定一致）
        theta = acosd(max(-1, min(1, R_rel(3, 3))));
        phi = atan2d(R_rel(2, 3), R_rel(1, 3));

        pose_seq(:, k) = [t_world; phi; theta];
    end
end

function [R, O] = build_world_frame(p1, p2, p3)
% 世界坐标系：三点均值为原点，p2指向p3为x轴，
% p1到x轴的垂线（由垂足指向p1）为y轴，z轴由右手定则确定
    O = (p1 + p2 + p3) / 3;
    x = (p3 - p2) / norm(p3 - p2);
    foot = p2 + x * dot(p1 - p2, x);  % p1 在x轴上的垂足
    y = (p1 - foot) / norm(p1 - foot);
    z = cross(x, y);
    R = [x y z];
end

function [R, O] = build_frame(p1, p2, p3)
% 测量子坐标系：原点为三点均值，x轴为 p3 指向 p2 的方向，
% y轴为 p1 到x轴的垂线方向（由p1指向垂足），z轴由右手定则确定
    O = (p1 + p2 + p3) / 3;
    x = (p2 - p3) / norm(p2 - p3);
    foot = p3 + x * dot(p1 - p3, x);  % p1 在x轴上的垂足
    y = (foot - p1) / norm(foot - p1);
    z = cross(x, y);
    R = [x y z];
end

function [ids, pts] = read_points_file(file_name)
% 读取测量点文件，返回点编号列向量 ids 与 3×n 坐标矩阵 pts
    fid = fopen(file_name, 'r');
    if fid < 0
        error('无法打开文件: %s', file_name);
    end
    ids = [];
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
        tok = regexp(strtrim(tokens{1}), '[a-zA-Z]+(\d+)', 'tokens', 'once');
        ids(end+1, 1) = str2double(tok{1}); %#ok<AGROW>
        pts(:, end+1) = str2double(strtrim(tokens(2:4))).'; %#ok<AGROW>
    end
    fclose(fid);
end
