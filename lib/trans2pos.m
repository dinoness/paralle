function pos = trans2pos(T)
%TRANS2POS 齐次变换矩阵 → 位姿向量 (pos2trans 的逆)
%   pos = [x; y; z; phi; theta]
%     phi   — z轴在xOy平面投影与x轴的夹角 (deg)
%     theta — z轴与世界z轴 [0;0;1] 的夹角 (deg)

t = T(1:3, 4);
z_axis = T(1:3, 3);

theta = atan2d(norm(z_axis(1:2)), z_axis(3));
phi = atan2d(z_axis(2), z_axis(1));

pos = [t; phi; theta];
end