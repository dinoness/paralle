% 参考文献：c02 并联机构的运动学误差建模及参数可辨识性分析_孔令雨
clear
addpath(genpath('./lib'));
%% 参数集
%--------parameter3--------
T = readtable('parameters.xlsx', 'Range', 'A2:B12');
paras = table2array(T(:, 2));
l_max = paras(1);
l_min = paras(2);  % 670
R1 = paras(3);  % 550
R2 = paras(4);  % 500
H = paras(5);  % 0
r1 = paras(6);  % 100
r2 = paras(7);  % 80
h = paras(8);  % 10

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

% -----end-parameter3------

% ----- input data ------
Pos_m_seq = [0.01;-0.01;-600.02;0.001;-0.001];  % line=5 colum=n
Pos_ref_seq = [0;30;-600;0;0];  % line=5 colum=n  角度的单位是° **一列为一组**
seq_len = length(Pos_ref_seq(1, :));
Pos_err_seq = zeros(5, seq_len);  % 位姿估计误差，优化的目标
Pos_delta_seq = zeros(5, seq_len);  % 位姿扰动序列
% ----- end input data ------


%% 运动学逆解部分（包括求解被动关节角度）
i_seq = 1;  % 循环变量，后期记得改=================
Pos_m = Pos_m_seq(:, i_seq);
Pos_ref = Pos_ref_seq(:, i_seq);
s_limb = zeros(3, 5);  % 支链的方向向量
p_pos = Pos_ref(1 : 3);
joint_para = zeros(6, 5);  % line = 关节数  colum = 支链数 存储格式：R-R-P-S，SPR支链为(Nan)-R-P-S
% ==========若零位下支链长度为l0，则输入的支链指令长度为l-l0==========


% 平台姿态
theta_plant = Pos_ref(5) / 180 * pi;  % 绕 theta
phi_plant = Pos_ref(4) / 180 * pi;  % 绕 phi
zb = [sin(theta_plant)*cos(phi_plant);
      sin(theta_plant)*sin(phi_plant);
      cos(theta_plant)];

ObB1 = B1 - p_pos;  
xb = cross(ObB1, zb) / norm(cross(ObB1, zb));
yb = cross(zb, xb);

R_plant = [xb yb zb];
for i = 1 : 5
    P_v(:, i) = R_plant * P_m(:, i);
    P(:, i) = P_v(:, i) + p_pos;
end

% 计算支链方向与长度
for i_limb = 1 : length(P_v(1, :))
    vAa = p_pos + P_v(:, i_limb) - B(:, i_limb);
    len_vAa = norm(vAa);
    s_limb(:, i_limb) = vAa / len_vAa;
    joint_para(3, i_limb) = len_vAa;
end

% ================================== 求解SPR支链 ==================================



% ================================== 求解UPS支链 ==================================
% ======================== 零位建模 ========================
% U关节轴线方向
o13 = zeros(1, 3);
joint_u_angle_tilt = 155 / 180 * pi;
% 固定坐标到关节1的转移矩阵
joint_1_z = [sin(joint_u_angle_tilt) * cos(limb_dir(2));
            sin(joint_u_angle_tilt) * sin(limb_dir(2));
            cos(joint_u_angle_tilt)];
joint_1_x = [cos(limb_dir(2) - pi/2);
            sin(limb_dir(2) - pi/2);
            0];
joint_1_y = cross(joint_1_z, joint_1_x);
R01 = [joint_1_x joint_1_y joint_1_z];
t01 = B2;
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
l0 = 600;
R34 = eye(3);
t34 = [0;0;l0];
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
p_x = [0;sin(pi-limb_dir(2));cos(pi-limb_dir(2))];
p_y = [0;sin(3/2*pi-limb_dir(2));cos(3/2*pi-limb_dir(2))];
R_p = [p_x p_y p_z];
t_p = [0;0;r1];
T_p = [R_p t_p;o13 1];

% 验证位姿矩阵用
% T01*T12*T23*T34*T45*T56*T_p
% T01*T12
T_e0 = T01*T12*T23*T34*T45*T56*T_p;  % 参考零位下的转移矩阵
T_ref = [R_plant p_pos; o13 1];

% 局部指数基公式
zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
zeta_p = [0;0;0;0;0;1];  % 平移基底

% ======================== 求解部分 ========================
% ============ 移动副 ============
q3 = joint_para(3, 2) - l0;
T_zeta3 = exp_se3(zeta_p*q3);
T_xi3 = (T01*T12*T23) * T_zeta3 / (T01*T12*T23);
% ============ 虎克铰 ============
% r_u和p_u，首先此处都是在基底坐标系{B}下的表示，r_u直接选择旋转前的数据（通过参考零位+横向位移获得），p_u直接选择旋转后的数据（旋转后的球铰链位置）
% 转轴的方向也是在零位下，基底坐标系{B}下的表示
% 至于为什么论文中的算法能成功目前还没想清楚
% **每个旋量坐标都被定义为零位状态下的值**
r_s = P(:, 2);
temp = T01*T12*T23*T34;
r_s0 = temp(1:3, 4);
r_u = B(:, 2);
p_u = [eye(3) zeros(3,1)] * T_xi3 * [r_s0;1];  % 旋转前，将零位的点变换到经过T_xi3的位置
q_u = [eye(3) zeros(3,1)] * T_ref / T_e0 * [r_s0;1]; % 旋转后，将零位的点逆变换到经过T_xi1,2,3的位置
% q_u = r_s; % 旋转后
u_u = p_u - r_u;
v_u = q_u - r_u;

axis1 = R01*[0;0;1];
axis2 = R01*R12*[0;0;1];

alpha = ((axis1'*axis2)*axis2'*u_u - axis1'*v_u)/((axis1'*axis2)^2-1);
beta = ((axis1'*axis2)*axis2'*v_u - axis2'*u_u)/((axis1'*axis2)^2-1);
gamma2 = (u_u'*u_u-alpha^2-beta^2-2*alpha*beta*axis1'*axis2)/((cross(axis1,axis2))'*(cross(axis1,axis2)));
gamma = sqrt(gamma2);

z_u = alpha*axis1 + beta*axis2 + gamma*cross(axis1, axis2);

v_u_ = v_u - axis1*axis1'*v_u;
u_u_ = u_u - axis2*axis2'*u_u;
z1_u_ = z_u - axis1*axis1'*z_u;
z2_u_ = z_u - axis2*axis2'*z_u;
q1 = atan(-(axis1'*cross(v_u_, z1_u_))/(v_u_'*z1_u_));
q2 = atan((axis2'*cross(u_u_, z2_u_))/(u_u_'*z2_u_));

T_zeta1 = exp_se3(zeta_r*q1);
T_zeta2 = exp_se3(zeta_r*q2);
T_xi1 = T01 * T_zeta1 / T01;
T_xi2 = (T01*T12) * T_zeta2 / (T01*T12);

% ============ 球铰 ============
% 在计算u_s和v_s中，论文中的r_s，用零位中的球关节代替
% 至于原因，还没想通，目前认为是在计算q_s中，将位姿逆变换回了最接近零位的状态
r_w0 = T_e0(1:3, 4);
p_s = r_w0;
q_s = [eye(3) zeros(3,1)] * ((T_xi1*T_xi2*T_xi3) \ T_ref / T_e0) * ([r_w0;1]);
u_s = p_s - r_s0;
v_s = q_s - r_s0;

axis4 = R01*R12*R23*R34*[0;0;1];
axis5 = R01*R12*R23*R34*R45*[0;0;1];

alpha = ((axis4'*axis5)*axis5'*u_s - axis4'*v_s)/((axis4'*axis5)^2-1);
beta = ((axis4'*axis5)*axis5'*v_s - axis5'*u_s)/((axis4'*axis5)^2-1);
gamma2 = (u_s'*u_s-alpha^2-beta^2-2*alpha*beta*axis4'*axis5)/(cross(axis4, axis5)'*cross(axis4, axis5));
gamma = -sqrt(gamma2);  % 正负号判断

z_s = alpha*axis4 + beta*axis5 + gamma*cross(axis4, axis5);

v_s_ = v_s - axis4*axis4'*v_s;
u_s_ = u_s - axis5*axis5'*u_s;
z1_s_ = z_s - axis4*axis4'*z_s;
z2_s_ = z_s - axis5*axis5'*z_s;
q4 = atan(-(axis4'*cross(v_s_, z1_s_))/(v_s_'*z1_s_));
q5 = atan((axis5'*cross(u_s_, z2_s_))/(u_s_'*z2_s_));

T_zeta4 = exp_se3(zeta_r*q4);
T_zeta5 = exp_se3(zeta_r*q5);
T_xi4 = (T01*T12*T23*T34) * T_zeta4 / (T01*T12*T23*T34);
T_xi5 = (T01*T12*T23*T34*T45) * T_zeta5 / (T01*T12*T23*T34*T45);



% t_rw = [10;0;r1];  % 在绝对坐标系中，动平台往下10
% T_rw = [R_p t_rw;o13 1];
% temp = T01*T12*T23*T34*T45*T56*T_rw;
% r_w_ = temp(1:3, 4);
r_w_ = p_pos + R_plant*[0;0;-10];
axis6 = R01*R12*R23*R34*R45*R56*[0;0;1];
u_w = r_w_ - r_s;
v_w = [eye(3) zeros(3,1)] * ((T_xi1*T_xi2*T_xi3*T_xi4*T_xi5) \ T_ref / T_e0) * ([r_w_;1] - [r_s;1]);
v_w_ = v_w - axis6 * axis6' * v_w;
u_w_ = u_w - axis6 * axis6' * u_w;
q6 = atan((axis6'*cross(u_w_,v_w_))/(u_w_'*v_w_));


T_zeta6 = exp_se3(zeta_r*q6);
T_xi6 = (T01*T12*T23*T34*T45*T56) * T_zeta6 / (T01*T12*T23*T34*T45*T56);

T_total1 = T01*T_zeta1*T12*T_zeta2*T23*T_zeta3*T34  % 验证虎克铰解算
r_s
T_total2 = T01*T_zeta1*T12*T_zeta2*T23*T_zeta3*T34*T_zeta4*T45*T_zeta5*T56*T_zeta6*T_p  % 验证球铰解算
T_ref

% 正确数据，q2 = -0.6435 or -36.87°
% T_total第三列为[0.5196;0.3;-0.8]
% q4 = 1.3642e-15
% q5 = -0.6435


% 先把球铰链最后一个关节变换写好，然后准备组会PPT，确定不开后，再接着写SPR关节解算


% T01*T_zeta1*T12*T_zeta2*T23*T_zeta3*T34*T45*T56*T_p
% T_xi1 * T_xi2 * T_xi3 * T_e0
