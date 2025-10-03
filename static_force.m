clear

% Parameters
R1 = 800;
R2 = 600;
H = 20;
r1 = 100;
r2 = 80;
h = 100;
Lmax = 1000;
Lmin = 650;

% Force
g = 9.8;
% F_ex = [0; 7*g; 0; 0; 0; 0];  % y
F_ex = [7*g; 0; 0; 0; 0; 7*g*50];  % x
% F_ex = [0; 0; 7*g; 0; 0; 0];  % z


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
pos_plant = [0; 0; -600];
alpha_plant = 0 / 180 * pi;
beta_plant = 0 / 180 * pi;
gamma_plant = 0 / 180 * pi;

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


% Length and Poster of limbs
s_limb = zeros(3, 5);
l_limb = zeros(1, 5);

for i = 1:5
    v_limb = pos_plant + R_plant * P_m(:, i) - B(:, i);
    l_limb(i) = norm(v_limb);
    s_limb(:, i) = v_limb / l_limb(i);
end

% Jacobian Matrix
% part 1
x_m = [1; 0; 0];
J1 = [s_limb (R_plant * x_m)];

% part 2
P_v = zeros(3, 5);  % 只变换了方向，没变换起点
P = zeros(3, 5);    % 末端点坐标
for i = 1 : 5
    P_v(:, i) = R_plant * P_m(:, i);
    P(:, i) = P_v(:, i) + pos_plant;
end
% cross常量取1表示列向量叉乘
% 这两个J2取哪个还不知道，复杂点的是通过约束方程得出的雅可比矩阵直接换算过来的，简单点的是通过受力分析得来的
% J2 = [cross(P_v, s_limb, 1)  ...
%      (cross(P_v(:, 1 ), (R_plant * x_m)) + l_limb(1)*cross((R_plant * x_m), s_limb(:, 1)))];
J2 = [cross(P_v, s_limb, 1)  ...
      cross(P_v(:, 1 ), (R_plant * x_m))];

J = [J1' J2'];

% force solve
F_in = (J')\F_ex;  % 注意雅克比矩阵需要转置，具体为何见推导过程



% plot
fig = figure('Color', [1 1 1]);
plot3(B(1, :), B(2, :), B(3, :), 'o', 'Color', '#FF7F50');
hold on
plot3(P(1, :), P(2, :), P(3, :), 'o', 'Color', '#32CD32');
for i = 1 : 5
    plot3([B(1, i) P(1, i)], [B(2, i) P(2, i)], [B(3, i) P(3, i)], '-', 'Color', '#4682B4');
end

axis equal

% draw force
% 画arrow需要横向量
K_force_draw = 10;

p_st = pos_plant;
p_ed = p_st + F_ex(1:3) * K_force_draw;
arrow3(p_st', p_ed', 'r');

for i = 1 : 5
    if abs(F_in(i)) > 0.01
        p_st = P(:, i);
        p_ed = p_st + s_limb(:, i) * F_in(i) * K_force_draw;
        arrow3(p_st', p_ed', 'r');
    end
end

if abs(F_in(6)) > 0.01
    p_st = P(:, 1);
    p_ed = p_st + R_plant * x_m * F_in(6) * K_force_draw;
    arrow3(p_st', p_ed', 'r');
end

grid on
xlabel("x");
ylabel("y")
zlabel("z")






