% 目的：通过螺旋定理验证机构自由度

clear

% Parameters
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

pos_plant = [0; 0; -555];  % 后面作图用，不参与空间搜索
alpha_plant = paras(9) / 180 * pi;  % 绕 x
beta_plant = paras(10) / 180 * pi;  % 绕 y
gamma_plant = paras(11) / 180 * pi;  % 绕 z

% Force
g = 9.8;
% F_ex = [0; 10*g; 0; 0; 0; 0];  % y
% F_ex = [10*g; 0; 0; 0; 0; 0];  % x
F_ex = [0; 0; -30*g; 0; 0; 0];  % z
% F_ex = [0; -3.5*g; -sqrt(3)*3.5*g; 0; 0; 0];


% Points
B1 = [R1*cos(pi/2);   R1*sin(pi/2);   0];
B2 = [R1*cos(7*pi/6); R1*sin(7*pi/6); 0];
B3 = [R1*cos(-pi/6);  R1*sin(-pi/6);  0];
B4 = [R2*cos(pi/6);   R2*sin(pi/6);   H];
B5 = [R2*cos(5*pi/6); R2*sin(5*pi/6); H];
B = [B1 B2 B3 B4 B5];

P1_m = [r1*cos(pi/2);   r1*sin(pi/2);   0];
P2_m = [r1*cos(7*pi/6); r1*sin(7*pi/6); 0];
P3_m = [r1*cos(-pi/6);  r1*sin(-pi/6);  0];
P4_m = [r2*cos(pi/6);   r2*sin(pi/6);   h];
P5_m = [r2*cos(5*pi/6); r2*sin(5*pi/6); h];
P_m = [P1_m P2_m P3_m P4_m P5_m];


% Position and Posture of Move Plant
% pos_plant = [0; 0; -600];
% alpha_plant = 0 / 180 * pi;
% beta_plant = 0 / 180 * pi;
% gamma_plant = 0 / 180 * pi;

Rx = [1                0                 0;
      0 cos(alpha_plant) -sin(alpha_plant);
      0 sin(alpha_plant)  cos(alpha_plant)];

Ry = [ cos(beta_plant) 0 sin(beta_plant);
                     0 1               0;
      -sin(beta_plant) 0 cos(beta_plant)];

Rz = [cos(gamma_plant) -sin(gamma_plant) 0;
      sin(gamma_plant)  cos(gamma_plant) 0;
                     0                0  1];

R_plant = Rz * Ry * Rx;

P_v = zeros(3, 5);  % 只变换了方向，没变换起点
P = zeros(3, 5);    % 末端点坐标
for i = 1 : 5
    P_v(:, i) = R_plant * P_m(:, i);
    P(:, i) = P_v(:, i) + pos_plant;
end


% Length and Poster of limbs
s_limb = zeros(3, 5);
l_limb = zeros(1, 5);

for i = 1:5
    v_limb = pos_plant + R_plant * P_m(:, i) - B(:, i);
    l_limb(i) = norm(v_limb);
    s_limb(:, i) = v_limb / l_limb(i);
end

% 1号支链
s11 = [1 0 0];
s12 = [0 1 0];
s13 = [0 0 1];
s14 = s_limb(:, 1)';
s15 = R_plant(:, 1)';

S11 = [s11 cross(B1', s11)];
S12 = [s12 cross(B1', s12)];
S13 = [s13 cross(B1', s13)];
S14 = [0 0 0 s14];
S15 = [s15 cross(P(:, 1)', s15)];

% 以球铰作为坐标原点
% S11 = [s11 0 0 0];
% S12 = [s12 0 0 0];
% S13 = [s13 0 0 0];
% S14 = [0 0 0 s14];
% S15 = [s15 cross(s14*l_limb(1), s15)];

TP1 = [S11; S12; S13; S14; S15];
TEMP1 = (null(TP1))';
C1 = [TEMP1(4:6) TEMP1(1:3)]
% 两种结果不一样？？？


% 2号支链
angle_tilt = 65 / 180 * pi;
s21 = [cos(pi/6)*cos(angle_tilt) sin(pi/6)*cos(angle_tilt) sin(angle_tilt)];
s23 = s_limb(:, 2)';
s22 = cross(s21, s23)/norm(cross(s21, s23));
s24 = [1 0 0];
s25 = [0 1 0];
s26 = [0 0 1];


S21 = [s21 cross(B2', s21)];
S22 = [s22 cross(B2', s22)];
S23 = [0 0 0 s23];
S24 = [s24 cross(P(:, 2)', s24)];
S25 = [s25 cross(P(:, 2)', s25)];
S26 = [s26 cross(P(:, 2)', s26)];

TP2 = [S21; S22; S23; S24; S25; S26];
