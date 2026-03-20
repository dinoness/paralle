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

% function S = skew(w)
% % skew    Returns the 3x3 skew‑symmetric matrix of a 3‑vector.
% %   S = skew(w) generates the matrix such that S * v = w × v for any v.
%     S = [0,    -w(3),  w(2);
%          w(3),  0,    -w(1);
%         -w(2),  w(1),  0];
% end

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
Te0 = T0*T1*Tp;

T1zeta = exp_se3(zeta_r * q1);
T2zeta = exp_se3(zeta_r * q2);

Tx1 = T0 * T1zeta / T0;
Tx2 = (T0*T1) * T2zeta / (T0*T1);

xi1 = log_se3(T0 * exp_se3(zeta_r) / T0);
xi2 = log_se3((T0*T1) * exp_se3(zeta_r) / (T0*T1));
Te = Tx1 * Tx2 * Te0;

J1 = xi1;
J2 = adjoint(Tx1,xi2);
disp("全局法")
disp([J1 J2])
% 全局法，旋量本身要零位的，而算伴随矩阵要带上现有位姿的

TJ1 = T0 * exp_se3(zeta_r) / T0;
j1 = log_se3(TJ1);


TJ2 = (T0*T1zeta*T1) * exp_se3(zeta_r) / (T0*T1zeta*T1);
j2 = log_se3(TJ2);

disp("局部法")
disp([j1 j2])

Jb1 = adjoint(trans_inv(Te), J1);
Jb2 = adjoint(trans_inv(Te), J2);

disp("body jacobian")
disp([Jb1 Jb2])

tol = 1e-5;
q1d = q1 + tol;
q2d = q2;
% 扰动正解
T1zetad = exp_se3(zeta_r * q1d);
T2zetad = exp_se3(zeta_r * q2d);
Tx1d = T0 * T1zetad / T0;
Tx2d = (T0*T1) * T2zetad / (T0*T1);
Ted = T0*T1zetad*T1*T2zetad*Tp;
% Ted = Tx1d * Tx2d * Te0;
% 求解扰动量
err = log_se3(Ted/Te) / tol
err2 = log_se3(Te\Ted) / tol


function j = adjoint(T, xi)
    R = T(1:3, 1:3);
    t = T(1:3, 4);
    j = [R zeros(3,3); skew(t)*R R] * xi;
end

function T_inv = trans_inv(T)
    R = T(1:3, 1:3);
    t = T(1:3, 4);
    T_inv = [R' -R'*t; 0 0 0 1];    
end

function S = skew(w)
% skew    Returns the 3x3 skew‑symmetric matrix of a 3‑vector.
%   S = skew(w) generates the matrix such that S * v = w × v for any v.
    S = [0,    -w(3),  w(2);
         w(3),  0,    -w(1);
        -w(2),  w(1),  0];
end


%% from deepseek
% % 参数
% l1 = 1;
% l2 = 0.8;
% theta0 = [0.5; 0.8];  % 测试位形
% delta = 1e-7;          % 扰动大小

% % 解析物体雅可比
% J_ana = jacobian_body(theta0, l1, l2);

% % 当前位姿
% T0 = fkine(theta0, l1, l2);

% % 初始化数值雅可比
% J_num = zeros(6,2);

% for j = 1:2
%     theta = theta0;
%     theta(j) = theta(j) + delta;
%     T = fkine(theta, l1, l2);
%     % 计算相对变换，取对数得物体旋量
%     dT = T / T0;   % 等价于 T * inv(T0)
%     dxi = log_se3(dT);
%     J_num(:,j) = dxi / delta;
% end

% % 比较
% err = norm(J_ana - J_num, 'fro') / norm(J_ana, 'fro');
% fprintf('Relative Frobenius error: %e\n', err);

% % 逐列比较
% for j = 1:2
%     col_err = norm(J_ana(:,j) - J_num(:,j)) / norm(J_ana(:,j));
%     fprintf('Column %d relative error: %e\n', j, col_err);
% end


% function Jb = jacobian_body(theta, l1, l2)
%     % 输入 theta: [θ1; θ2]
%     % 输出 Jb: 6x2 物体雅可比

%     xi1 = [0;0;1; 0;0;0];
%     xi2 = [0;0;1; 0;l1;0];
%     M = [1 0 0 l1+l2; 0 1 0 0; 0 0 1 0; 0 0 0 1];

%     % 计算关节2的空间旋量
%     T1 = exp_se3(xi1 * theta(1));
%     Ad_T1 = adjoint(T1);
%     xi2_space = Ad_T1 * xi2;

%     % 空间雅可比
%     Js = [xi1, xi2_space];

%     % 末端位姿
%     T = fkine(theta, l1, l2);
%     Ad_Tinv = adjoint(invT(T));  % Ad_{T^{-1}}
%     Jb = Ad_Tinv * Js;
% end

% function Ad = adjoint(T)
%     R = T(1:3,1:3);
%     t = T(1:3,4);
%     Ad = [R, zeros(3,3); skew(t)*R, R];
% end

% function Tinv = invT(T)
%     R = T(1:3,1:3);
%     t = T(1:3,4);
%     Tinv = [R', -R'*t; 0 0 0 1];
% end

% function T = fkine(theta, l1, l2)
%     % theta: [θ1; θ2] 关节角度（弧度）
%     % 返回末端位姿矩阵 4x4

%     xi1 = [0;0;1; 0;0;0];
%     xi2 = [0;0;1; 0;l1;0];
%     M = [1 0 0 l1+l2; 0 1 0 0; 0 0 1 0; 0 0 0 1];

%     T = exp_se3(xi1 * theta(1)) * exp_se3(xi2 * theta(2)) * M;
% end




% function T = exp_se3(xi)
%     % xi = [ω; v] 6x1 旋量坐标
%     w = xi(1:3);
%     v = xi(4:6);
%     theta = norm(w);
%     if theta < 1e-10
%         R = eye(3);
%         p = v;
%     else
%         w_hat = skew(w/theta);  % 单位方向向量的反对称矩阵
%         R = eye(3) + sin(theta)*w_hat + (1-cos(theta))*w_hat^2;
%         % 平移部分： (Iθ + (1-cosθ)w_hat + (θ-sinθ)w_hat^2) * v/θ
%         p = (eye(3)*theta + (1-cos(theta))*w_hat + (theta-sin(theta))*w_hat^2) * v / theta;
%     end
%     T = [R, p; 0 0 0 1];
% end

% function S = skew(w)
%     S = [0, -w(3), w(2);
%          w(3), 0, -w(1);
%          -w(2), w(1), 0];
% end