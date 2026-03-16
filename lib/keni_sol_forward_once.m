function Te_seq = keni_sol_forward_once(joint_q, limb_dir, B, r1, r2, l0_seq, P_m)
% 正向运动学，求解一次
% input： 关节结构参数，关节位姿
% output：各直连平台位姿

Te_seq = zeros(4,4,5);



joint_u_angle_tilt = 155 / 180 * pi;
o13 = zeros(1, 3);
% 局部指数基公式
zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
zeta_p = [0;0;0;0;0;1];  % 平移基底




% SPR支链
% ======================== 零位建模 ========================
q1 = joint_q(:, 1);
q15 = q1(2);
q14 = q1(3);
q11 = q1(4);
q12 = q1(5);
q13 = q1(6);

% 固定坐标到关节1的转移矩阵
j1_1z = [1;0;0];
j1_1x = [0;1;0];
j1_1y = [0;0;1];
R1_1 = [j1_1x j1_1y j1_1z];
t1_1 = B(:,1);
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

T1_zeta1 = exp_se3(zeta_r*q11);
T1_zeta2 = exp_se3(zeta_r*q12);
T1_zeta3 = exp_se3(zeta_r*q13);
T1_zeta4 = exp_se3(zeta_p*q14);
T1_zeta5 = exp_se3(zeta_r*q15);

Te_seq(:,:,1) = T1_1*T1_zeta1*T1_2*T1_zeta2*T1_3*T1_zeta3*T1_4*T1_zeta4*T1_5*T1_zeta5*T1_p;


% UPS支链
for i_limb = 2 : 5
    if i_limb <= 3
        r = r1;
    else
        r = r2;
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

    q = joint_q(:, i_limb);
    T_zeta1 = exp_se3(zeta_r*q(1));
    T_zeta2 = exp_se3(zeta_r*q(2));
    T_zeta3 = exp_se3(zeta_p*q(3));
    T_zeta4 = exp_se3(zeta_r*q(4));
    T_zeta5 = exp_se3(zeta_r*q(5));
    T_zeta6 = exp_se3(zeta_r*q(6));

    Te_seq(:,:,i_limb) = T01*T_zeta1*T12*T_zeta2*T23*T_zeta3*T34*T_zeta4*T45*T_zeta5*T56*T_zeta6*T_p;
end



end