% 单层检索并作图
clear
addpath(genpath('./lib'));
fprintf('>>>= start (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));

%% 参数集
%--------parameter3--------
unit_para = 0.001;  % 0.001表示m，1表示mm

T = readtable('parameters.xlsx', 'Range', 'A2:B13');
paras = table2array(T(:, 2));
l_max = paras(1)*unit_para;
l_min = paras(2)*unit_para;  % 670
R1 = paras(3)*unit_para;  % 550
R2 = paras(4)*unit_para;  % 500
H = paras(5)*unit_para;  % 0
r1 = paras(6)*unit_para;  % 100
r2 = paras(7)*unit_para;  % 80
h = paras(8)*unit_para;  % 10
% L_tool = paras(12)*unit_para;
L_tool = 0;  % ===========================================

% limb_dir = [pi/2; 7*pi/6; -pi/6; pi/6; 5*pi/6];
limb_dir = [pi/2; 7*pi/6; -pi/6; pi/2-deg2rad(45); pi/2+deg2rad(45)];
% limb_dir = [pi/2; deg2rad(-90-45); deg2rad(-90+45); deg2rad(45); deg2rad(135)];
B1 = [R1*cos(limb_dir(1)); R1*sin(limb_dir(1)); 0];
B2 = [R1*cos(limb_dir(2)); R1*sin(limb_dir(2)); 0];
B3 = [R1*cos(limb_dir(3)); R1*sin(limb_dir(3)); 0];
B4 = [R2*cos(limb_dir(4)); R2*sin(limb_dir(4)); H];
B5 = [R2*cos(limb_dir(5)); R2*sin(limb_dir(5)); H];
B = [B1 B2 B3 B4 B5];

% move plant parameter
limb_dir_move = limb_dir;
P1_m = [r1*cos(limb_dir_move(1)); r1*sin(limb_dir_move(1)); L_tool];
P2_m = [r1*cos(limb_dir_move(2)); r1*sin(limb_dir_move(2)); L_tool];
P3_m = [r1*cos(limb_dir_move(3)); r1*sin(limb_dir_move(3)); L_tool];
P4_m = [r2*cos(limb_dir_move(4)); r2*sin(limb_dir_move(4)); L_tool + h];
P5_m = [r2*cos(limb_dir_move(5)); r2*sin(limb_dir_move(5)); L_tool + h];
P_m = [P1_m P2_m P3_m P4_m P5_m];
P_v = zeros(3, 5);  % 只变换了方向，没变换起点
P = zeros(3, 5);    % 末端点坐标

% ball screw dir vector
ball_screw_dir_angle_deg = [145 145 145  145 145;   % 与Z轴夹角
                            -90  30 150 -150 -30];  % 与x轴夹角
ball_screw_dir_angle = ball_screw_dir_angle_deg / 180 * pi;
ball_vector = zeros(3, 5);
ball_vector_world = zeros(3, 5);

static_joint_dir_angle_deg = [ 0  0  0  0  0;
                              -90  30 150 -150 -30];
static_joint_dir_angle = static_joint_dir_angle_deg / 180 * pi;
static_joint_vector = zeros(3, 5);

for i_ball = 1 : 5
    ball_vector(1, i_ball) = sin(ball_screw_dir_angle(1, i_ball)) * cos(ball_screw_dir_angle(2, i_ball));
    ball_vector(2, i_ball) = sin(ball_screw_dir_angle(1, i_ball)) * sin(ball_screw_dir_angle(2, i_ball));
    ball_vector(3, i_ball) = cos(ball_screw_dir_angle(1, i_ball));
    

    static_joint_vector(1, i_ball) = sin(static_joint_dir_angle(1, i_ball)) * cos(static_joint_dir_angle(2, i_ball));
    static_joint_vector(2, i_ball) = sin(static_joint_dir_angle(1, i_ball)) * sin(static_joint_dir_angle(2, i_ball));
    static_joint_vector(3, i_ball) = cos(static_joint_dir_angle(1, i_ball));
end

% 局部指数基公式
zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
zeta_p = [0;0;0;0;0;1];  % 平移基底
Pos_ref_seq = [0*unit_para;0*unit_para;-600*unit_para;0;0];  % 角度单位是°
T_ref = pos2trans(Pos_ref_seq(:, 1), B);
l0 = 600;
l0_seq = [l0;l0;l0;l0;l0];
joint_u_angle_tilt = 155 / 180 * pi;
% -----end-parameter3------

[p_seq, xi_seq] = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);
Omega = [zeros(3,3) eye(3); eye(3) zeros(3,3)];
tol = 1e-9;


seq_x = -200 : 10 : 200;
seq_y = -200 : 10 : 200;
seq_z = -1000 : 20 : -400;
seq_phi = -30 : 10 : 30;
seq_theta = -20 : 5 : 20;
seq_x = seq_x * unit_para;
seq_y = seq_y * unit_para;
seq_z = seq_z * unit_para;
l_seq_x = length(seq_x);
l_seq_y = length(seq_y);
l_seq_z = length(seq_z);
l_seq_phi = length(seq_phi);
l_seq_theta = length(seq_theta);

data_plot = zeros(3, l_seq_x*l_seq_y);

num_count = 0;
sum_TF = 0;

parfor ix = 1 : l_seq_x
    for iy = 1 : l_seq_y
        for iz = 1 : l_seq_z
            for iphi = 1 : l_seq_phi
                for itheta = 1 : l_seq_theta
                    % 先判断是否在工作区间中
                    pos_flag = 0;  % 位置可达标志位
                    Pos_ref = [seq_x(ix); seq_y(iy); seq_z(iz); seq_phi(iphi); seq_theta(itheta)];
                    T_ref = pos2trans(Pos_ref, B);  

                    % 当前位形下平台关键点坐标
                    R_plant = T_ref(1:3, 1:3);
                    pos_plant = T_ref(1:3, 4);
                    P_v = R_plant * P_m;
                    P = P_v + pos_plant;

                    ball_vector_world = R_plant * ball_vector;

                    for j = 1 : length(P_v(1, :))
                        vAa = pos_plant + P_v(:, j) - B(:, j);
                        len_vAa = norm(vAa);
                        s_limb = vAa / len_vAa;
                        
                        if (len_vAa >= l_min)&&(len_vAa <= l_max)  % ===========支链长度条件===========
                            angle_limb_scew = acos(dot(s_limb, ball_vector_world(:, j)));  % 支链与球铰轴线夹角
                            angle_limb_scew_deg = angle_limb_scew / pi * 180;
                            if(angle_limb_scew_deg <= 30 )  % ===========关节角度条件===========
                                pos_flag = pos_flag + 1;
                            end
                        end
                    end



                    if pos_flag < length(P_v(1, :))
                        continue;
                    end

                    joint_q0 = keni_sol_inverse(T_ref, B, l0_seq, P_m, p_seq);

                    % 变量定义
                    U = zeros(6, 6, 5);     % 每条支链的关节运动旋量（不足6列的部分补零）
                    ST = zeros(6, 5);       % transmission wrench screw，按列存放 [F; M]
                    SI = zeros(6, 5);       % 输入旋量（主动P副）
                    SO = zeros(6, 5);       % output twist screw，按列存放 [omega; v]
                    SC = zeros(6, 1);       % 机构公共约束力旋量
                    passive_recip_dim = zeros(5, 1);

                    %% 1) 计算各支链关节运动旋量与 TWS
                    for i_limb = 1 : 5
                        [Ui, active_idx, passive_idx] = spr4ups_build_limb_twists(p_seq, joint_q0, i_limb, zeta_r, zeta_p);
                        U(:, 1:size(Ui, 2), i_limb) = Ui;

                        % 主动关节输入旋量（P副）
                        SI(:, i_limb) = screw_unitize(Ui(:, active_idx), 'twist');

                        % TWS：对主动P副，传递力旋量可直接取为沿支链轴线、并通过平台连接点Pi的纯力
                        limb_dir_now = Ui(4:6, active_idx);
                        ST(:, i_limb) = build_line_wrench(P(:, i_limb), limb_dir_now);
                        ST(:, i_limb) = screw_unitize(ST(:, i_limb), 'wrench');

                        % 调整符号，使 TWS 与输入旋量做功为正
                        if SI(:, i_limb).' * Omega * ST(:, i_limb) < 0
                            ST(:, i_limb) = -ST(:, i_limb);
                        end

                        % 记录“仅对被动关节取倒易”的维数，便于检查SPR与UPS支链是否符合理论预期
                        passive_recip_dim(i_limb) = size(svd_nullspace(Ui(:, passive_idx).' * Omega, tol), 2);
                    end

                    %% 2) 计算公共约束力旋量 SC（来源于 SPR 支链）
                    U1 = U(:, 1:5, 1);
                    SC_basis = svd_nullspace(U1.' * Omega, tol);
                    if isempty(SC_basis)
                        error('SPR 支链对全部关节取倒易后未得到公共约束力旋量，请检查机构参数或是否处于奇异位形。');
                    end
                    SC = screw_unitize(SC_basis(:, 1), 'wrench');

                    % 统一SC的符号：若其力部与平台x轴反向，则翻转，便于结果解释
                    if norm(SC(1:3)) > tol
                        if dot(SC(1:3), R_plant(:, 1)) < 0
                            SC = -SC;
                        end
                    end

                    %% 3) 计算各支链 OTS
                    % 定义：当其余4个主动P副全部锁定时，平台在公共约束 SC + 其余4条支链 TWS 的作用下
                    % 剩余的唯一允许瞬时运动，即为该支链对应的 OTS。
                    for i_limb = 1 : 5
                        other_idx = setdiff(1:5, i_limb);
                        W_locked = [SC, ST(:, other_idx)];
                        SO_basis = svd_nullspace(W_locked.' * Omega, tol);

                        if isempty(SO_basis)
                            error('第 %d 条支链的 OTS 未求得，请检查是否处于奇异位形。', i_limb);
                        end

                        SO(:, i_limb) = screw_unitize(SO_basis(:, 1), 'twist');

                        % 统一OTS符号，使其与本支链TWS的瞬时功率为正
                        if ST(:, i_limb).' * Omega * SO(:, i_limb) < 0
                            SO(:, i_limb) = -SO(:, i_limb);
                        end
                    end

                    %% 4) 倒易性检查
                    err_TWS_passive = zeros(5, 1);
                    err_SC = norm(U1.' * Omega * SC);
                    err_OTS = zeros(5, 1);
                    err_TWS_active_power = zeros(5, 1);

                    for i_limb = 1 : 5
                        if i_limb == 1
                            passive_idx = [1 2 3 5];
                            Ui = U(:, 1:5, i_limb);
                        else
                            passive_idx = [1 2 4 5 6];
                            Ui = U(:, 1:6, i_limb);
                        end

                        err_TWS_passive(i_limb) = norm(Ui(:, passive_idx).' * Omega * ST(:, i_limb));
                        err_TWS_active_power(i_limb) = SI(:, i_limb).' * Omega * ST(:, i_limb);

                        other_idx = setdiff(1:5, i_limb);
                        W_locked = [SC, ST(:, other_idx)];
                        err_OTS(i_limb) = norm(W_locked.' * Omega * SO(:, i_limb));
                    end

                    %% 4.5) 计算式(3)-(4)中的 MF / TF

                    MF_num = zeros(5,1);   % 分子 |ST o SI|
                    MF_den = zeros(5,1);   % 分母 |ST o SI|_max
                    MF = zeros(5,1);       % 归一化 MF

                    TF_num = zeros(5,1);   % 分子 |ST o SO|
                    TF_den = zeros(5,1);   % 分母 |ST o SO|_max
                    TF = zeros(5,1);       % 归一化 TF

                    for i_limb = 1 : 5
                        % reciprocal product
                        MF_num(i_limb) = abs(ST(:, i_limb).' * Omega * SI(:, i_limb));
                        TF_num(i_limb) = abs(ST(:, i_limb).' * Omega * SO(:, i_limb));

                        % -------- MF denominator --------
                        % 当前机构中 SI 是主动P副的纯平移单位 twist，
                        % ST 是沿支链轴线的单位 pure wrench，
                        % 二者最大 reciprocal product = 1
                        MF_den(i_limb) = 1.0;

                        % -------- TF denominator --------
                        % 对单位 wrench screw 和单位 twist screw，
                        % 最大 reciprocal product 取 sqrt((h_w + h_t)^2 + d^2)
                        % 这里 ST 是 pure wrench，因此 h_w = 0
                        TF_den(i_limb) = max_reciprocal_wrench_twist(ST(:, i_limb), SO(:, i_limb), tol);

                        % 归一化
                        MF(i_limb) = MF_num(i_limb) / MF_den(i_limb);

                        if TF_den(i_limb) > tol
                            TF(i_limb) = TF_num(i_limb) / TF_den(i_limb);
                        else
                            TF(i_limb) = 0;
                        end
                    end

                    num_count = num_count + 1;
                    sum_TF = sum_TF + min(TF);


                end
            end
        end
    end
end
fprintf('>>>= done (%s) , average TF = %.6f =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')), sum_TF/num_count);




%% Functions
function vmax = max_reciprocal_wrench_twist(w, t, tol)
% 计算 |w o t|_max
% w: 6x1 wrench [F; M]
% t: 6x1 twist  [omega; v]

    [sw, rw, hw, is_couple] = wrench_axis_pitch(w, tol);
    [st, rt, ht, is_translation] = twist_axis_pitch(t, tol);

    if is_couple
        error('当前函数暂不处理纯偶力 wrench 作为 TWS 的情况。');
    end

    % 若 OTS 是纯平移，则单位力与单位平移的最大 reciprocal product 为 1
    if is_translation
        vmax = 1.0;
        return;
    end

    % 两 screw 轴线间的最短距离 d
    d = line_distance(rw, sw, rt, st, tol);

    % 最大 reciprocal product
    vmax = sqrt((hw + ht)^2 + d^2);
end


function [s, r, h, is_couple] = wrench_axis_pitch(w, tol)
% 从 wrench [F; M] 提取轴线方向 s、轴上一点 r、pitch h
% 约定：M = r x F + h F

    F = w(1:3);
    M = w(4:6);

    if norm(F) < tol
        is_couple = true;
        s = M / norm(M);
        r = [0;0;0];
        h = inf;
        return;
    end

    is_couple = false;

    FF = dot(F, F);
    s = F / norm(F);
    h = dot(F, M) / FF;

    % 取轴线上到原点的垂足点
    r = cross(F, M) / FF;
end


function [s, r, h, is_translation] = twist_axis_pitch(t, tol)
% 从 twist [omega; v] 提取轴线方向 s、轴上一点 r、pitch h
% 对应 twist 形式：v = -omega x r + h * omega

    omega = t(1:3);
    v = t(4:6);

    if norm(omega) < tol
        is_translation = true;
        s = v / norm(v);
        r = [0;0;0];
        h = inf;
        return;
    end

    is_translation = false;

    ww = dot(omega, omega);
    s = omega / norm(omega);
    h = dot(omega, v) / ww;

    % 取轴线上到原点的垂足点
    r = cross(omega, v) / ww;
end


function d = line_distance(r1, s1, r2, s2, tol)
% 两空间直线的最短距离
% line1: r = r1 + lambda*s1
% line2: r = r2 + mu*s2

    c = cross(s1, s2);
    nc = norm(c);

    if nc < tol
        % 平行线
        d = norm(cross(s1, (r2 - r1)));
    else
        % 斜交线
        d = abs(dot((r2 - r1), c)) / nc;
    end
end