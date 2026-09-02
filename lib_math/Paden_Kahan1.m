function q = Paden_Kahan1(u,v,axis)
% 求解Paden_Kahan第一个子问题，解算两轴旋转的逆运动
% 输入：
%   u：起始向量
%   v：终止向量
%   axis1：轴线
% 输出：
%   q：转角（弧度）

v_ = v - axis * axis' * v;
u_ = u - axis * axis' * u;
q = atan((axis'*cross(u_,v_))/(u_'*v_));

end