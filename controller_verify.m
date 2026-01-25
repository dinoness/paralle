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


pos_plant = [0; 0; -700];
theta = 10 / 180 * pi;
phi = 0 / 180 * pi;


% static plant
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

zb = [cos(phi)*sin(theta); sin(phi)*sin(theta); cos(theta)];
ObB1 = - pos_plant + B1;
xb = cross(ObB1, zb) / norm(cross(ObB1, zb));
yb = cross(zb, xb);

R_plant = [xb yb zb];
P_v = zeros(3, 5);  % 只变换了方向，没变换起点
for i = 1 : 5
    P_v(:, i) = R_plant * P_m(:, i);
    
end

l_limb = zeros(5, 1);
for i = 1 : 5
    l_limb(i) = norm(pos_plant + P_v(:, i) - B(:, i));
end

l_limb
