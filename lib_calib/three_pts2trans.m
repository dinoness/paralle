function T = three_pts2trans(p1, p2, p3)
%THREE_PTS_TO_TRANS 由三个不共线点计算坐标系M到空间原点的转移矩阵
%   输入：
%     p1, p2, p3 — 3×1 向量，三个点的三维坐标 [x; y; z]
%   输出：
%     T — 4×4 齐次变换矩阵，满足 p_world = T * p_M
%
%   坐标系M定义：
%     原点OM：三点重心
%     z轴：三点所在平面的法向量，与世界z轴 [0;0;1] 夹角 ≤ 90°
%     y轴：点1指向点2和点3的中点
%     x轴：右手定则 x = y × z

    OM = (p1 + p2 + p3) / 3;
    mid_23 = (p2 + p3) / 2;

    y = mid_23 - p1;
    y = y / norm(y);

    z = cross(p2 - p1, p3 - p1);
    z = z / norm(z);
    if z(3) < 0
        z = -z;
    end

    x = cross(y, z);
    x = x / norm(x);
    y = cross(z, x);

    R = [x, y, z];
    T = [R, OM; 0, 0, 0, 1];
end
