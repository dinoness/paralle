function joint_q = keni_sol_inverse(Te, B, l0_seq, P_m, p_seq)
% 运动学逆解
% input：结构参数，平台位姿
% output：各直连关节量

joint_q = zeros(6,5);


o13 = zeros(1, 3);
% 局部指数基公式
zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
zeta_p = [0;0;0;0;0;1];  % 平移基底



P = zeros(3, 5);
p_pos = Te(1:3, 4);
R_plant = Te(1:3, 1:3);
limb_length = zeros(5, 1);
% 计算支链方向与长度
for i_limb = 1 : 5
    vAa = p_pos + R_plant * P_m(:, i_limb) - B(:, i_limb);
    P(:, i_limb) = R_plant * P_m(:, i_limb) + p_pos;
    len_vAa = norm(vAa);
    limb_length(i_limb) = len_vAa;
end


% SPR支链
T1 = zeros(4,4,6);
for i_joint = 1 : 6
    T1(:,:,i_joint) = exp_se3(p_seq(:, i_joint));
end
T1_1 = T1(:,:,1);
T1_2 = T1(:,:,2);
T1_3 = T1(:,:,3);
T1_4 = T1(:,:,4);
T1_5 = T1(:,:,5);
T1_p = T1(:,:,6);

T1e0 = T1_1*T1_2*T1_3*T1_4*T1_5*T1_p;

% ======================== 求解部分 ========================
% ============ 移动副 ============
q14 = limb_length(1) - l0_seq(1);
T1_zeta4 = exp_se3(zeta_p*q14);
T1_xi4 = (T1_1*T1_2*T1_3*T1_4) * T1_zeta4 / (T1_1*T1_2*T1_3*T1_4);

% ============ 球副 ============
% 前两轴
temp = T1_1*T1_2*T1_3*T1_4*T1_5;
r1r0 = temp(1:3, 4);
r1s = B(:, 1);
p1s = [eye(3) zeros(3,1)] * T1_xi4 * [r1r0;1];
q1s = [eye(3) zeros(3,1)] * (Te / T1e0) * [r1r0;1];

axis11 = Rotate_axis(T1_1);
axis12 = Rotate_axis(T1_1*T1_2);

u1s = p1s - r1s;
v1s = q1s - r1s;
[q11,q12] = Paden_Kahan2(u1s,v1s,axis11,axis12,1);

T1_zeta1 = exp_se3(zeta_r*q11);
T1_zeta2 = exp_se3(zeta_r*q12);
T1_xi1 = T1_1 * T1_zeta1 / T1_1;
T1_xi2 = (T1_1*T1_2) * T1_zeta2 / (T1_1*T1_2);

% 第三轴
dtemp = T1_1*T1_2*T1_3*T1_4*T1_5*[eye(3) [0;0;10];o13 1];
r1w = dtemp(1:3,4);
axis13 = Rotate_axis(T1_1*T1_2*T1_3);

p1w = [eye(3) zeros(3,1)] * T1_xi4 * ([r1w;1]);
q1w = [eye(3) zeros(3,1)] * ((T1_xi1*T1_xi2) \ Te / T1e0) * ([r1w;1]);
u1w = p1w - r1s;
v1w = q1w - r1s;
q13 = Paden_Kahan1(u1w, v1w, axis13);
T1_zeta3 = exp_se3(zeta_r*q13);
T1_xi3 = (T1_1*T1_2*T1_3) * T1_zeta3 / (T1_1*T1_2*T1_3);


% %验证用
% cal_r = T1_1*T1_zeta1*T1_2*T1_zeta2*T1_3*T1_zeta3*T1_4*T1_zeta4*T1_5
% real_r = P(:, 1)


% ============ 回转副 ============
axis15 = Rotate_axis(T1_1*T1_2*T1_3*T1_4*T1_5);
temp = T1_1*T1_2*T1_3*T1_4*T1_5*T1_p;
r1p = temp(1:3, 4);
p1p = r1p;
q1p = [eye(3) zeros(3,1)] * ((T1_xi1*T1_xi2*T1_xi3*T1_xi4) \ Te / T1e0) * [r1p;1];
u1p = p1p - r1r0;
v1p = q1p - r1r0;
q15 = Paden_Kahan1(u1p,v1p,axis15);

joint_q(:, 1) = [q11;q12;q13;q14;q15;0];

% UPS
for i_limb = 2 : 5

    T2 = zeros(4, 4, 7);
    for i_joint = 1 : 7
        T2(:,:,i_joint) = exp_se3(p_seq(:, (i_limb-1)*7 + (i_joint-1)));
    end
    T01 = T2(:,:,1);
    T12 = T2(:,:,2);
    T23 = T2(:,:,3);
    T34 = T2(:,:,4);
    T45 = T2(:,:,5);
    T56 = T2(:,:,6);
    T_p = T2(:,:,7);

    T_e0 = T01*T12*T23*T34*T45*T56*T_p;  % 参考零位下的转移矩阵

    % ============ 移动副 ============
    q3 = limb_length(i_limb) - l0_seq(i_limb);
    T_zeta3 = exp_se3(zeta_p*q3);
    T_xi3 = (T01*T12*T23) * T_zeta3 / (T01*T12*T23);

    % ============ 虎克铰 ============
    r_s = P(:, i_limb);
    temp = T01*T12*T23*T34;
    r_s0 = temp(1:3, 4);
    r_u = B(:, i_limb);
    p_u = [eye(3) zeros(3,1)] * T_xi3 * [r_s0;1];  % 旋转前，将零位的点变换到经过T_xi3的位置
    q_u = [eye(3) zeros(3,1)] * Te / T_e0 * [r_s0;1]; % 旋转后，将零位的点逆变换到经过T_xi1,2,3的位置
    % q_u = r_s; % 旋转后
    u_u = p_u - r_u;
    v_u = q_u - r_u;

    axis1 = Rotate_axis(T01);
    axis2 = Rotate_axis(T01*T12);

    [q1,q2] = Paden_Kahan2(u_u, v_u, axis1, axis2, 1);

    T_zeta1 = exp_se3(zeta_r*q1);
    T_zeta2 = exp_se3(zeta_r*q2);
    T_xi1 = T01 * T_zeta1 / T01;
    T_xi2 = (T01*T12) * T_zeta2 / (T01*T12);

    % ============ 球铰 ============
    r_w0 = T_e0(1:3, 4) + [0;0;P_m(3,i_limb)];
    p_s = r_w0;
    q_s = [eye(3) zeros(3,1)] * ((T_xi1*T_xi2*T_xi3) \ Te / T_e0) * ([r_w0;1]);
    u_s = p_s - r_s0;
    v_s = q_s - r_s0;

    axis4 = Rotate_axis(T01*T12*T23*T34);
    axis5 = Rotate_axis(T01*T12*T23*T34*T45);

    [q4,q5] = Paden_Kahan2(u_s, v_s, axis4, axis5, -1);

    T_zeta4 = exp_se3(zeta_r*q4);
    T_zeta5 = exp_se3(zeta_r*q5);
    T_xi4 = (T01*T12*T23*T34) * T_zeta4 / (T01*T12*T23*T34);
    T_xi5 = (T01*T12*T23*T34*T45) * T_zeta5 / (T01*T12*T23*T34*T45);

    r_w_ = p_pos + R_plant*[0;0;-10];
    axis6 = Rotate_axis(T01*T12*T23*T34*T45*T56);
    u_w = r_w_ - r_s;
    v_w = [eye(3) zeros(3,1)] * ((T_xi1*T_xi2*T_xi3*T_xi4*T_xi5) \ Te / T_e0) * ([r_w_;1] - [r_s;1]);
    q6 = Paden_Kahan1(u_w,v_w,axis6);

    joint_q(: ,i_limb) = [q1;q2;q3;q4;q5;q6];
end


end

function axis = Rotate_axis(T)
    axis = T(1:3,3);
end