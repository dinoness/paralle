function [pose_seq, T_seq] = calib_pts2pose_seq(data_dir)
%CALIB_PTS2POSE_SEQ 处理标定测量数据，生成动平台表面坐标系位姿序列
%   输入：
%     data_dir — 数据目录，包含 t1.txt / t2.txt / t3.txt 三个文件，
%                每行格式：行号(制表符)点编号, x, y, z（默认 'calib_p'）
%   输出：
%     pose_seq — 5×n 位姿序列，每列 [x; y; z; phi; theta]，
%                为动平台表面坐标系在世界坐标系下的位姿，
%                长度单位 mm，角度单位 deg
%     T_seq    — 4×4×n 齐次变换矩阵序列，由测量点构建的坐标系直接
%                生成，为动平台表面坐标系在世界坐标系下的表示，
%                平移单位 mm
%
%   处理流程：
%     1. 每个 txt 文件的首行与末行单独取出取平均，作为该文件的参考点，
%        其余行为测量点；
%     2. 参考坐标系：原点为三个参考点的均值，x轴为 t3参考点 指向 t2参考点
%        的方向，y轴为原点指向x轴的垂线方向，z轴由右手定则确定；
%     3. 相同编号的三个测量点以同样方法构建测量子坐标系，求其在参考
%        坐标系下的表示；
%     4. 世界坐标系：坐标轴与参考坐标系平行，参考坐标系原点在世界系中
%        的位置为 (0, 0, -837.8) mm；
%     5. 动平台表面坐标系：坐标轴与测量子坐标系平行，其原点在测量子坐标
%        系中的位置为 (0, 0, 37.8) mm；最终输出动平台表面坐标系在世界
%        坐标系下的位姿序列。

    if nargin < 1 || isempty(data_dir)
        data_dir = 'calib_p';
    end

    O_ref_in_world = [0; 0; -837.8];  % 参考坐标系原点在世界系中的位置 (mm)
    O_surf_in_meas = [0; 0; 37.8];    % 动平台表面原点在测量子坐标系中的位置 (mm)

    [ids1, pts1] = read_points_file(fullfile(data_dir, 't1.txt'));
    [ids2, pts2] = read_points_file(fullfile(data_dir, 't2.txt'));
    [ids3, pts3] = read_points_file(fullfile(data_dir, 't3.txt'));

    % 1. 首末行取平均作为参考点，其余为测量点
    ref1 = mean(pts1(:, [1 end]), 2);
    ref2 = mean(pts2(:, [1 end]), 2);
    ref3 = mean(pts3(:, [1 end]), 2);

    ids1 = ids1(2:end-1);  pts1 = pts1(:, 2:end-1);
    ids2 = ids2(2:end-1);  pts2 = pts2(:, 2:end-1);
    ids3 = ids3(2:end-1);  pts3 = pts3(:, 2:end-1);

    % 2. 参考坐标系
    [R_ref, O_ref] = build_frame(ref1, ref2, ref3);

    % 3. 相同编号的测量点构建测量子坐标系
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

        % 测量子坐标系在参考坐标系下的表示
        R_rel = R_ref.' * R_sub;
        t_rel = R_ref.' * (O_sub - O_ref);

        % 4-5. 测量子坐标系 -> 动平台表面坐标系 -> 世界坐标系
        % T_world_surf = T_world_ref * T_ref_meas * T_meas_surf
        % 其中 T_world_ref 与 T_meas_surf 均为纯平移（坐标轴相互平行）
        t_world = O_ref_in_world + t_rel + R_rel * O_surf_in_meas;
        T_seq(:, :, k) = [R_rel, t_world; 0 0 0 1];

        % 姿态不变，由相对旋转矩阵的 z 轴列恢复 phi/theta
        % （与 pos2trans 约定一致）
        theta = acosd(max(-1, min(1, R_rel(3, 3))));
        phi = atan2d(R_rel(2, 3), R_rel(1, 3));

        pose_seq(:, k) = [t_world; phi; theta];
    end
end

function [R, O] = build_frame(p1, p2, p3)
% 三点构建坐标系：原点为三点均值，x轴为 p3 指向 p2 的方向，
% y轴为原点指向x轴的垂线方向，z轴由右手定则确定
    O = (p1 + p2 + p3) / 3;
    x = (p2 - p3) / norm(p2 - p3);
    foot = p3 + x * dot(O - p3, x);  % 原点在x轴上的垂足
    y = (foot - O) / norm(foot - O);
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
