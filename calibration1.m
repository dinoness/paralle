% 输入：一系列测量位姿Pm，理论位姿P0，理想结构参数
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

alpha_plant = paras(9) / 180 * pi;  % 绕 x
beta_plant = paras(10) / 180 * pi;  % 绕 y
gamma_plant = paras(11) / 180 * pi;  % 绕 z

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
Pos_m_seq = [];  % line=5 colum=n
Pos_ref_seq = [];  % line=5 colum=n
seq_len = length(Pos_ref_seq(1, :));
Pos_err_seq = zeors(5, seq_len);  % 位姿估计误差，优化的目标
Pos_delta_seq = zeors(5, seq_len);  % 位姿扰动序列
% ----- end input data ------

% ----- output data ------
d_B = zeros(3, 5);
d_P = zeros(3, 5);

% ----- end output data ------


for i_seq = 1 : seq_len
    Pos_m = Pos_m_seq(:, i_seq);
    Pos_ref = Pos_ref_seq(:, i_seq);

    s_limb = zeros(3, 5);  % 支链的方向向量

    % 平台姿态
    alpha_plant = Pos_ref(4);
    beta_plant = Pos_ref(5);
    gamma_plant = Pos_ref(6);
    Rx = [1                  0                   0;
          0   cos(alpha_plant)   -sin(alpha_plant);
          0   sin(alpha_plant)    cos(alpha_plant)];

    Ry = [ cos(beta_plant)   0   sin(beta_plant);
                         0   1                 0;
          -sin(beta_plant)   0   cos(beta_plant)];

    Rz = [cos(gamma_plant) -sin(gamma_plant)   0;
          sin(gamma_plant)  cos(gamma_plant)   0;
                         0                0    1];

    R_plant = Rz * Ry * Rx;
    for i = 1 : 5
        P_v(:, i) = R_plant * P_m(:, i);
        P(:, i) = P_v(:, i) + pos_plant;
    end

    % 计算支链方向
    for i_limb = 1 : length(P_v(1, :))
        vAa = vt + P_v(:, i_limb) - B(:, i_limb);
        len_vAa = norm(vAa);
        s_limb(:, i_limb) = vAa / len_vAa;
    end

    % 求解末端扰动
    Pos_delta = zeros(5, 1);
    % J_1  line = 5, colum = 6
    J_1 = [s_limb(1)' cross(R_plant*B(:, 1),s_limb(1))';
           s_limb(2)' cross(R_plant*B(:, 2),s_limb(2))';
           s_limb(3)' cross(R_plant*B(:, 3),s_limb(3))';
           s_limb(4)' cross(R_plant*B(:, 4),s_limb(4))';
           s_limb(5)' cross(R_plant*B(:, 5),s_limb(5))'];
    % J_p  line = 6, colum = 6


    b_sol = zeros(5, 1);
    for i_limb = 1 : length(b_sol)
        b_sol(i_limb) = -s_limb(:, i_limb)' * R_plant * d_P(:, i_limb) ...
                        + s_limb(:, i_limb)' * d_B(:, i_limb);
    end

    Pos_delta = pinv(J_1) * b_sol;  % 用伪逆来解方程 不一定可行
    Pos_delta_seq(:, i_seq) = Pos_delta;

end

% 后续的问题：如果用这个来做了，怎么做误差的最小化呢？
% 换句话说：求解算法需要的格式是什么？求解非线性最小二乘？

