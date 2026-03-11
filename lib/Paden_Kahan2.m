function [q1, q2] = Paden_Kahan2(u,v,axis1,axis2,sign)
% 求解Paden_Kahan第二个子问题，解算两轴旋转的逆运动
% 输入：
%   p：起始向量
%   q：终止向量
%   r：轴线的交点
%   axis1：第一轴轴线
%   axis2：第二轴轴线
%   sign_：用于gamma开方的正负号，之后可以找方法进行自动判断，待优化
% 输出：
%   q1：第一轴的转角（弧度）
%   q2：第二轴的转角（弧度）

if sign > 0
    sign_ = 1;
else
    sign_ = -1;
end

alpha = ((axis1'*axis2)*axis2'*u - axis1'*v)/((axis1'*axis2)^2-1);
beta = ((axis1'*axis2)*axis2'*v - axis2'*u)/((axis1'*axis2)^2-1);
gamma2 = (u'*u-alpha^2-beta^2-2*alpha*beta*axis1'*axis2)/((cross(axis1,axis2))'*(cross(axis1,axis2)));
gamma = sign_*sqrt(gamma2);

z = alpha*axis1 + beta*axis2 + gamma*(cross(axis1,axis2));

v_ = v - axis1*axis1'*v;
u_ = u - axis2*axis2'*u;
z1_ = z - axis1*axis1'*z;
z2_ = z - axis2*axis2'*z;

q1 = atan(-(axis1'*cross(v_, z1_))/(v_'*z1_));
q2 = atan((axis2'*cross(u_, z2_))/(u_'*z2_));

end