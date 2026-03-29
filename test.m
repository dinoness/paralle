% === 测试脚本 ===
clear;
addpath(genpath('./AI_code'));

%--------struct parameter--------
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
B_ = [B1 B2 B3 B4 B5];

% move plant parameter
P1_m = [r1*cos(pi/2);   r1*sin(pi/2);   0];
P2_m = [r1*cos(7*pi/6); r1*sin(7*pi/6); 0];
P3_m = [r1*cos(-pi/6);  r1*sin(-pi/6);  0];
P4_m = [r2*cos(pi/6);   r2*sin(pi/6);   h];
P5_m = [r2*cos(5*pi/6); r2*sin(5*pi/6); h];
P_m = [P1_m P2_m P3_m P4_m P5_m];

params.B = B_;
params.A_local = P_m;

% 1. 定义机构几何参数 (根据你的CAD模型替换这些虚拟值)
% 基座铰链点 B (3x5) - 假设分布在圆上
% Rb = 0.5; 
% angles = linspace(0, 2*pi, 6); angles(end) = [];
% params.B = [Rb*cos(angles); Rb*sin(angles); zeros(1,5)];

% 动平台局部铰链点 A_local (3x5)
% Ra = 0.3;
% params.A_local = [Ra*cos(angles); Ra*sin(angles); zeros(1,5)];

% 2. 定义当前位姿 (例如：Home Position)
% 绕X轴旋转10度，Z轴平移0.8米
theta = deg2rad(0);
R_test = [1, 0, 0; 0, cos(theta), -sin(theta); 0, sin(theta), cos(theta)];
P_test = [0; 0; -0.6];
T_mat = eye(4);
T_mat(1:3, 1:3) = R_test;
T_mat(1:3, 4) = P_test;

% 3. 定义外力 W_ext (施加在动平台中心)
% 例如：受到Z方向 -100N 的负载，以及X方向 5Nm 的扭矩
W_ext = [0; 0; -294; 0; 0; 1]; 

% 4. 设定SPR约束条件
% 假设SPR支链在位置1，其R关节在动平台上，局部轴线为[1;0;0]
R_axis_local = [1; 0; 0];
params.R_axis_global = R_test * R_axis_local; % 转换到全局系

% 5. 调用核心函数求解
[tau, fc, G] = SPR_4UPS_Statics(T_mat, W_ext, params);

% 6. 输出结果
disp('5个P关节的主动驱动力 (N):');
disp(tau);
disp('SPR支链提供的约束反力矩 (Nm):');
disp(fc);