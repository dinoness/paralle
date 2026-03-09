% 参考文献：c02 并联机构的运动学误差建模及参数可辨识性分析_孔令雨


addpath(genpath('./lib'));
%% 参数集
%--------parameter3--------
T = readtable('parameters.xlsx', 'Range', 'A2:B12');
paras = table2array(T(:, 2));
l_max = paras(1);
l_min = paras(2);  % 670
R1 = paras(3);  % 800
R2 = paras(4);  % 600
H = paras(5);  % 20
r1 = paras(6);  % 100
r2 = paras(7);  % 80
h = paras(8);  % 100


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
Pos_ref_seq = [0;0;-600;0;0];  % line=5 colum=n  角度的单位是° **一列为一组**
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
phi_plant = paras(10) / 180 * pi;  % 绕 phi
zb = [cos(theta_plant)*sin(phi_plant);
        cos(theta_plant)*cos(phi_plant);
        sin(theta_plant)];

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
    joint_para(3, i) = len_vAa;
end

% 求解SPR支链



% 求解UPS支链
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


% 关节1到关节2的转移矩阵
joint_2_z = [1;0;0];
joint_2_x = [0;1;0];
joint_2_y = [0;0;1];
R12 = [joint_2_x joint_2_y joint_2_z];
t12 = [0;0;0];
T12 = [R12 t12;o13 1];


% 关节2到关节3的转移矩阵
l0 = 600;
joint_3_z = [1;0;0];
joint_3_x = [0;1;0];
joint_3_y = [0;0;1];
R23 = [joint_3_x joint_3_y joint_3_z];
t23 = [l0;0;0];
T23 = [R23 t23;o13 1];

% 关节3到关节4的转移矩阵
R34 = eye(3);
t34 = zeros(3,1);
T34 = [R34 t34;o13 1];

% 关节4到关节5的转移矩阵
joint_5_z = [0;1;0];
joint_5_x = [1;0;0];
joint_5_y = [0;0;-1];
R45 = [joint_5_x joint_5_y joint_5_z];
t45 = zeros(3,1);
T45 = [R45 t45;o13 1];

% 关节5到关节6的转移矩阵
joint_6_z = [1;0;0];
joint_6_x = [0;-1;0];
joint_6_y = [0;0;-1];
R56 = [joint_6_x joint_6_y joint_6_z];
t56 = zeros(3,1);
T56 = [R56 t56;o13 1];

% ============ 求解部分 ============
% 移动副
d_length = joint_para(3, 2) - l0;

% 虎克铰
r_s = P(:, 2);












