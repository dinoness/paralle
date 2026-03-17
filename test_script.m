addpath(genpath('./lib'));

%% 测试论文中的位姿数据
% p1 = [0;0;2.0944;0.2056;0.3560;0];
% p2 = [-pi/2;0;0;0;0;0];
% p3 = [0;0;0;0;0.5;0];
% p4 = [pi/2;0;0;-0.15;0;0];
% p5 = [0;pi/2;0;0;0;0];
% p6 = [1.5835;-0.9142;-1.5835;-0.0325;-0.1343;-0.1];
% p = [p1 p2 p3 p4 p5 p6];

% T = [1 0 0 0;
%      0 1 0 0;
%      0 0 1 0;
%      0 0 0 1];

% for i = 1 : 6
%     T0 = exp_se3(p(:,i))
%     T = T * T0;
% end

% fprintf("T = \n")
% disp(T)

%% 用建议虎克铰测试求解方法
% o13 = zeros(1, 3);
% % 固定坐标到关节1的转移矩阵
% joint_1_z = [1;0;0];
% joint_1_x = [0;1;0];
% joint_1_y = [0;0;1];
% R01 = [joint_1_x joint_1_y joint_1_z];
% t01 = [0;0;0];
% T01 = [R01 t01;o13 1];

% % 关节1到关节2的转移矩阵(绕z)
% joint_2_z = [1;0;0];
% joint_2_x = [0;0;-1];
% joint_2_y = [0;1;0];
% R12 = [joint_2_x joint_2_y joint_2_z];
% t12 = [0;0;0];
% T12 = [R12 t12;o13 1];

% % 关节2到动平台中心转移矩阵(绕z)
% p_z = [0;1;0];
% p_x = [-1;0;0];
% p_y = [0;0;1];
% R_p = [p_x p_y p_z];
% t_p = [0;5;0];
% T_p = [R_p t_p;o13 1];


% % 局部指数基公式
% zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
% zeta_p = [0;0;0;0;0;1];  % 平移基底


% % ======================
% % r = [0;0;0];
% % p = [0;0;5];  % 旋转前
% % phi = pi/3;
% % theta = pi/6;
% % q = 50*[sin(theta)*cos(phi);
% %          sin(theta)*sin(phi);
% %          cos(theta)];
% % u = p - r;
% % v = q - r;

% % axis1 = [1;0;0];
% % axis2 = [0;1;0];
% % axis1 = R01*[0;0;1];
% % axis2 = R01*R12*[0;0;1];
% % ======================
% r = [0;0;0];
% phi = pi/3;
% theta = pi/6;
% p = [sqrt(100^2-40^2-60^2);
%      40;
%      60];
% q = 50*[cos(pi/6);
%         sin(pi/6);
%         0];
% u = p - r;
% v = q - r;
% axis1 = [0;0;-1];
% axis2 = [0.5;cos(pi/6);0];

% alpha = ((axis1'*axis2)*axis2'*u - axis1'*v)/((axis1'*axis2)^2-1);
% beta = ((axis1'*axis2)*axis2'*v - axis2'*u)/((axis1'*axis2)^2-1);
% gamma2 = (u'*u-alpha^2-beta^2-2*alpha*beta*axis1'*axis2)/((skew(axis1)*axis2)'*(skew(axis1)*axis2));
% gamma = -sqrt(gamma2);

% z = alpha*axis1 + beta*axis2 + gamma*(skew(axis1)*axis2);


% v_ = v - axis1*axis1'*v;
% u_ = u - axis2*axis2'*u;
% z1_ = z - axis1*axis1'*z;
% z2_ = z - axis2*axis2'*z;
% q1 = atan(-(axis1'*(skew(v_)*z1_))/(v_'*z1_));
% q2 = atan((axis2'*(skew(u_)*z2_))/(u_'*z2_));

% T = T01*exp_se3(zeta_r*q1)*T12*exp_se3(zeta_r*q2)*T_p;

function S = skew(w)
% skew    Returns the 3x3 skew‑symmetric matrix of a 3‑vector.
%   S = skew(w) generates the matrix such that S * v = w × v for any v.
    S = [0,    -w(3),  w(2);
         w(3),  0,    -w(1);
        -w(2),  w(1),  0];
end

%% 测试雅克比矩阵求法
zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
zeta_p = [0;0;0;0;0;1];  % 平移基底

L1 = 10;
L2 = 5;
q1 = pi/6;
q2 = pi/4;

T0 = [eye(4)];
T1 = [eye(3) [L1;0;0]; [0 0 0 1]];
Tp = [eye(3) [L2;0;0]; [0 0 0 1]];

T1zeta = exp_se3(zeta_r * q1);
T2zeta = exp_se3(zeta_r * q2);

Tx1 = T0 * T1zeta / T0;
Tx2 = (T0*T1) * T2zeta / (T0*T1);

xi1 = log_se3(T0 * exp_se3(zeta_r) / T0);
xi2 = log_se3((T0*T1) * exp_se3(zeta_r) / (T0*T1));


J1 = xi1;
J2 = adjoint(Tx1,xi2);
disp("全局法")
disp([J1 J2])
% 全局法，旋量本身要零位的，而算伴随矩阵要带上现有位姿的

J1 = T0 * exp_se3(zeta_r) / T0;
j1 = log_se3(J1);
adjoint(T0, zeta_r);


J2 = (T0*T1zeta*T1) * exp_se3(zeta_r) / (T0*T1zeta*T1);
j2 = log_se3(J2);
adjoint((T0*T1zeta*T1), zeta_r);

disp("局部法")
disp([j1 j2])


function j = adjoint(T, xi)
    R = T(1:3, 1:3);
    t = T(1:3, 4);
    j = [R zeros(3,3); skew(t)*R R] * xi;
end