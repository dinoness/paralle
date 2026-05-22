function [p_seq, xi_seq] = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt)
% 给定几何参数，生成旋量序列
p_seq = zeros(6, 34); % 1-6,7-13,14-20,21-27,28-34
xi_seq = zeros(6, 34);  % 关节零位全局坐标
o13 = zeros(1, 3);

% 局部指数基公式
zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
zeta_p = [0;0;0;0;0;1];  % 平移基底


% SPR支链
% 固定坐标到关节1的转移矩阵
j1_1z = [1;0;0];
j1_1x = [0;1;0];
j1_1y = [0;0;1];
R1_1 = [j1_1x j1_1y j1_1z];
t1_1 = B(:, 1);
T1_1 = [R1_1 t1_1;o13 1];

% 关节1到关节2的转移矩阵
j1_2z = [-1;0;0];
j1_2x = [0;0;1];
j1_2y = [0;1;0];
R1_2 = [j1_2x j1_2y j1_2z];
t1_2 = zeros(3,1);
T1_2 = [R1_2 t1_2;o13 1];

% 关节2到关节3的转移矩阵
j1_3z = [0;-1;0];
j1_3x = [1;0;0];
j1_3y = [0;0;1];
R1_3 = [j1_3x j1_3y j1_3z];
t1_3 = zeros(3,1);
T1_3 = [R1_3 t1_3;o13 1];

% 关节3到关节4的转移矩阵
R1_4 = eye(3);
t1_4 = zeros(3,1);
T1_4 = [R1_4 t1_4;o13 1];

% 关节4到关节5的转移矩阵
j1_5z = [-1;0;0];
j1_5x = [0;0;1];
j1_5y = [0;1;0];
R1_5 = [j1_5x j1_5y j1_5z];
t1_5 = [0;0;l0_seq(1)];
T1_5 = [R1_5 t1_5;o13 1];

% 关节5到平台的转移矩阵
j1_pz = [-1;0;0];
j1_px = [0;0;-1];
j1_py = [0;-1;0];
R1_p = [j1_px j1_py j1_pz];
t1_p = [0;r1;0];
T1_p = [R1_p t1_p;o13 1];

p_seq(:, 1) = log_se3(T1_1);
p_seq(:, 2) = log_se3(T1_2);
p_seq(:, 3) = log_se3(T1_3);
p_seq(:, 4) = log_se3(T1_4);
p_seq(:, 5) = log_se3(T1_5);
p_seq(:, 6) = log_se3(T1_p);

xi_seq(:, 1) = adjoint_m(T1_1) * zeta_r;
xi_seq(:, 2) = adjoint_m(T1_1*T1_2) * zeta_r;
xi_seq(:, 3) = adjoint_m(T1_1*T1_2*T1_3) * zeta_r;
xi_seq(:, 4) = adjoint_m(T1_1*T1_2*T1_3*T1_4) * zeta_p;
xi_seq(:, 5) = adjoint_m(T1_1*T1_2*T1_3*T1_4*T1_5) * zeta_r;
% xi_seq(:, 6) = adjoint_m(T1_1*T1_2*T1_3*T1_4*T1_5*T1_p);

% UPS
for i_limb = 2 : 5
    if i_limb <= 3
        r = r1;
        % r_w0 = T_e0(1:3, 4);
    else
        r = r2;
        % r_w0 = T_e0(1:3, 4) + [0;0;P_m(5,i_limb)];
    end

    % 固定坐标到关节1的转移矩阵
    joint_1_z = [sin(joint_u_angle_tilt) * cos(limb_dir(i_limb));
                sin(joint_u_angle_tilt) * sin(limb_dir(i_limb));
                cos(joint_u_angle_tilt)];
    joint_1_x = [cos(limb_dir(i_limb) - pi/2);
                sin(limb_dir(i_limb) - pi/2);
                0];
    joint_1_y = cross(joint_1_z, joint_1_x);
    R01 = [joint_1_x joint_1_y joint_1_z];
    t01 = B(:, i_limb);
    T01 = [R01 t01;o13 1];


    % 关节1到关节2的转移矩阵(绕z)
    joint_2_z = [1;0;0];
    joint_2_x = [0;sin(pi-joint_u_angle_tilt);cos(pi-joint_u_angle_tilt)];
    joint_2_y = cross(joint_2_z, joint_2_x);
    R12 = [joint_2_x joint_2_y joint_2_z];
    t12 = [0;0;0];
    T12 = [R12 t12;o13 1];


    % 关节2到关节3的转移矩阵(绕z)
    joint_3_z = [1;0;0];
    joint_3_x = [0;0;1];
    joint_3_y = [0;-1;0];
    R23 = [joint_3_x joint_3_y joint_3_z];
    t23 = [0;0;0];
    T23 = [R23 t23;o13 1];

    % 关节3到关节4的转移矩阵(平移z)
    R34 = eye(3);
    t34 = [0;0;l0_seq(i_limb)];
    T34 = [R34 t34;o13 1];

    % 关节4到关节5的转移矩阵(绕z)
    joint_5_z = [-1;0;0];
    joint_5_x = [0;0;1];
    joint_5_y = [0;1;0];
    R45 = [joint_5_x joint_5_y joint_5_z];
    t45 = zeros(3,1);
    T45 = [R45 t45;o13 1];

    % 关节5到关节6的转移矩阵(绕z)
    joint_6_z = [0;1;0];
    joint_6_x = [1;0;0];
    joint_6_y = [0;0;-1];
    R56 = [joint_6_x joint_6_y joint_6_z];
    t56 = zeros(3,1);
    T56 = [R56 t56;o13 1];

    % 关节6到动平台中心转移矩阵(绕z)
    p_z = [-1;0;0];
    p_x = [0;sin(pi-limb_dir(i_limb));cos(pi-limb_dir(i_limb))];
    p_y = [0;sin(3/2*pi-limb_dir(i_limb));cos(3/2*pi-limb_dir(i_limb))];
    R_p = [p_x p_y p_z];
    t_p = [P_m(3,i_limb);0;r];
    T_p = [R_p t_p;o13 1];

    p_seq(:, 7*(i_limb-1) + 0) = log_se3(T01);
    p_seq(:, 7*(i_limb-1) + 1) = log_se3(T12);
    p_seq(:, 7*(i_limb-1) + 2) = log_se3(T23);
    p_seq(:, 7*(i_limb-1) + 3) = log_se3(T34);
    p_seq(:, 7*(i_limb-1) + 4) = log_se3(T45);
    p_seq(:, 7*(i_limb-1) + 5) = log_se3(T56);
    p_seq(:, 7*(i_limb-1) + 6) = log_se3(T_p);

    xi_seq(:, 7*(i_limb-1) + 0) = adjoint_m(T01) * zeta_r;
    xi_seq(:, 7*(i_limb-1) + 1) = adjoint_m(T01*T12) * zeta_r;
    xi_seq(:, 7*(i_limb-1) + 2) = adjoint_m(T01*T12*T23) * zeta_p;
    xi_seq(:, 7*(i_limb-1) + 3) = adjoint_m(T01*T12*T23*T34) * zeta_r;
    xi_seq(:, 7*(i_limb-1) + 4) = adjoint_m(T01*T12*T23*T34*T45) * zeta_r;
    xi_seq(:, 7*(i_limb-1) + 5) = adjoint_m(T01*T12*T23*T34*T45*T56) * zeta_r;
    % xi_seq(:, 7*(i_limb-1) + 6) = adjoint_m(T_p);

end


end