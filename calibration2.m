% 输入：一系列测量位姿Pm(五个变量)，理论位姿P0，理想结构参数
% 输出：结构参数的偏差值


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
Pos_ref_seq = [0;0;-600;0;0];  % line=5 colum=n  角度的单位是°
seq_len = length(Pos_ref_seq(1, :));
Pos_err_seq = zeros(5, seq_len);  % 位姿估计误差，优化的目标
Pos_delta_seq = zeros(5, seq_len);  % 位姿扰动序列
% ----- end input data ------

% ----- output data ------
d_B = zeros(3, 5);
d_P = zeros(3, 5);
d_L = zeros(1, 5);
for i = 1 : 3
    for j = 1 : 5
        d_P(i, j) = 0.02;
        d_B(i, j) = 0.03;
    end
end
% ----- end output data ------

% B1，B等参数之后可能会在迭代优化中发生改变======================================
e1 = [1; 0; 0];
e2 = [0; 1; 0];
e3 = [0; 0; 1];


for i_seq = 1 : seq_len
    Pos_m = Pos_m_seq(:, i_seq);
    Pos_ref = Pos_ref_seq(:, i_seq);

    s_limb = zeros(3, 5);  % 支链的方向向量
    p_pos = Pos_ref(1 : 3);

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

    % 计算支链方向
    for i_limb = 1 : length(P_v(1, :))
        vAa = p_pos + P_v(:, i_limb) - B(:, i_limb);
        len_vAa = norm(vAa);
        s_limb(:, i_limb) = vAa / len_vAa;
    end

    % 求解末端扰动
    Pos_delta = zeros(6, 1);
    % J1  line = 1, colum = 6
    J1 = [s_limb(:, 1)'  cross(R_plant*B(:, 1), s_limb(:, 1))'];
    J2 = [s_limb(:, 2)'  cross(R_plant*B(:, 2), s_limb(:, 2))'];
    J3 = [s_limb(:, 3)'  cross(R_plant*B(:, 3), s_limb(:, 3))'];
    J4 = [s_limb(:, 4)'  cross(R_plant*B(:, 4), s_limb(:, 4))'];
    J5 = [s_limb(:, 5)'  cross(R_plant*B(:, 5), s_limb(:, 5))'];
    Jc = [         -xb'  cross(             xb,         ObB1)'];
    % Jp  line = 6, colum = 6
    Jp = [J1;J2;J3;J4;J5;Jc];


    b_sol = zeros(6, 1);
    % 求解格式[d_Bi; d_Pi; d_Li]
    d_par1 = [d_B(:, 1); d_P(:, 1); d_L(1)];
    d_par2 = [d_B(:, 2); d_P(:, 2); d_L(2)];
    d_par3 = [d_B(:, 3); d_P(:, 3); d_L(3)];
    d_par4 = [d_B(:, 4); d_P(:, 4); d_L(4)];
    d_par5 = [d_B(:, 5); d_P(:, 5); d_L(5)];

    J11 = [s_limb(:, 1)'  -s_limb(:, 2)'*R_plant  1];
    J22 = [s_limb(:, 2)'  -s_limb(:, 2)'*R_plant  1];
    J33 = [s_limb(:, 3)'  -s_limb(:, 3)'*R_plant  1];
    J44 = [s_limb(:, 4)'  -s_limb(:, 4)'*R_plant  1];
    J55 = [s_limb(:, 5)'  -s_limb(:, 5)'*R_plant  1];
    Jcc = [         -xb'                     e1'  0];

    b_sol(1) = J11 * d_par1;
    b_sol(2) = J22 * d_par2;
    b_sol(3) = J33 * d_par3;
    b_sol(4) = J44 * d_par4;
    b_sol(5) = J55 * d_par5;
    b_sol(6) = Jcc * d_par1;
    

    Pos_delta = Jp \ b_sol;  % 求解位姿偏移量
    % Pos_delta_seq(:, i_seq) = Pos_delta;  % 角度表示的转换，求出的为alpha,beta,gamma形式欧拉角，而用于计算的是theta,phi倾角

end

% 当前有一个误差模型了(模型正确性待验证)
% 接下来的问题是如何将这个模型用于最小二乘优化

