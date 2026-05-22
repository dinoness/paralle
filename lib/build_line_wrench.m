function w = build_line_wrench(point_on_line, direction)
%BUILD_LINE_WRENCH 由“过一点、沿某方向”的直线构造纯力旋量 [F; M]
% 输入：
%   point_on_line : 3x1，作用线经过的一点
%   direction     : 3x1，作用线方向
% 输出：
%   w             : 6x1 wrench，格式为 [F; M]

    dir_norm = norm(direction);
    if dir_norm < 1e-12
        error('build_line_wrench: direction 的范数过小。');
    end

    e = direction / dir_norm;
    m = cross(point_on_line, e);
    w = [e; m];
end
