clear
addpath(genpath('./lib'));
%% 参数集
%--------parameter3--------
unit_para = 0.001;  % 0.001表示m，1表示mm

T = readtable('parameters.xlsx', 'Range', 'A2:B12');
paras = table2array(T(:, 2));
l_max = paras(1)*unit_para;
l_min = paras(2)*unit_para;  % 670
R1 = paras(3)*unit_para;  % 550
R2 = paras(4)*unit_para;  % 500
H = paras(5)*unit_para;  % 0
r1 = paras(6)*unit_para;  % 100
r2 = paras(7)*unit_para;  % 80
h = paras(8)*unit_para;  % 10

limb_dir = [pi/2; 7*pi/6; -pi/6; pi/6; 5*pi/6];
B1 = [R1*cos(pi/2);   R1*sin(pi/2);   0];
B2 = [R1*cos(7*pi/6); R1*sin(7*pi/6); 0];
B3 = [R1*cos(-pi/6);  R1*sin(-pi/6);  0];
B4 = [R2*cos(pi/6);   R2*sin(pi/6);   H];
B5 = [R2*cos(5*pi/6); R2*sin(5*pi/6); H];
B = [B1 B2 B3 B4 B5];

% move plant parameter
P1_m = [r1*cos(pi/2);   r1*sin(pi/2);   0];
P2_m = [r1*cos(7*pi/6); r1*sin(7*pi/6); 0];
P3_m = [r1*cos(-pi/6);  r1*sin(-pi/6);  0];
P4_m = [r2*cos(pi/6);   r2*sin(pi/6);   h];
P5_m = [r2*cos(5*pi/6); r2*sin(5*pi/6); h];
P_m = [P1_m P2_m P3_m P4_m P5_m];
P_v = zeros(3, 5);  % 只变换了方向，没变换起点
P = zeros(3, 5);    % 末端点坐标


% 局部指数基公式
zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
zeta_p = [0;0;0;0;0;1];  % 平移基底
Pos_ref_seq = [0*unit_para;0*unit_para;-600*unit_para;0;0];  % line=5 colum=n  角度的单位是° **一列为一组**
T_ref = pos2trans(Pos_ref_seq(:, 1), B);
l0 = 600;
l0_seq = [l0;l0;l0;l0;l0];
joint_u_angle_tilt = 155 / 180 * pi;

% -----end-parameter3------

p_seq = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

joint_q0 = keni_sol_inverse(T_ref, B, l0_seq, P_m, p_seq);

U = zeros(6, 6, 5);
ST = zeros(6,5);  % 传递力旋量
% SPR
U1 = zeros(6, 6);  % 运动副旋量

T1 = zeros(4,4,6);
for i_joint = 1 : 6
    T1(:,:,i_joint) = exp_se3(p_seq(:, i_joint));
end

T_temp = eye(4);
T_zeta = eye(4);
for i_joint = 1 : 5
    T_temp = T_temp * T_zeta * T1(:,:,i_joint);

    if i_joint ~= 4
        T_zeta = exp_se3(zeta_r * joint_q0(i_joint,1));
    else
        T_zeta = exp_se3(zeta_p * joint_q0(i_joint,1));
    end

    if i_joint ~= 4
        U1(:, i_joint) = adjoint_m(T_temp) * zeta_r;
    else
        U1(:, i_joint) = adjoint_m(T_temp) * zeta_p;
    end
end
U(:, :, 1) = U1;
Omega = [zeros(3,3) eye(3); eye(3) zeros(3,3)];
C1 = null(U1(:, 1:5)' * Omega);
ST(:, 1) = [U1(4:6,4);cross((B1 - Pos_ref_seq(1:3)),U1(4:6,4))];  % 纯力，将运动方向转移到上层，并以基座标为原点，叉乘上过支链点的矢量

% UPS
U2 = zeros(6, 6);
T2 = zeros(4,4,7);
for i_limb = 2 : 5
    for i_joint = 1 : 7
        T2(:,:,i_joint) = exp_se3(p_seq(:, (7*(i_limb-1) - 1) + i_joint));
    end

    for i_joint = 1 : 6
        T_temp = eye(4);
        for i_T = 1 : i_joint
            T_temp = T_temp * T_zeta * T2(:,:,i_T);
        end

        if i_joint ~= 3
            T_zeta = exp_se3(zeta_r * joint_q0(i_joint,i_limb));
        else
            T_zeta = exp_se3(zeta_p * joint_q0(i_joint,i_limb));
        end

        if i_joint ~= 3
            U2(:, i_joint) = adjoint_m(T_temp) * zeta_r;
        else
            U2(:, i_joint) = adjoint_m(T_temp) * zeta_p;
        end
    end

    U(:, :, i_limb) = U2;
    ST(:, i_limb) = [U2(4:6,3);cross((B(:,i_limb)- Pos_ref_seq(1:3)),U2(4:6,3))];
    % ST(:, i_limb) = null([U2(:, 1:2) U2(:, 4:6)]' * Omega);

end

% 目前的困难：无法验证驱动力旋量ST，以及约束旋量C1是否计算正确