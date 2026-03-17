function Te_seq = keni_sol_forward_once(joint_q, p_seq)
% 正向运动学，求解一次
% input： 关节结构参数，关节位姿
% output：各直连平台位姿

Te_seq = zeros(4,4,5);


% 局部指数基公式
zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
zeta_p = [0;0;0;0;0;1];  % 平移基底



% SPR支链
q1 = joint_q(:, 1);
q11 = q1(1);
q12 = q1(2);
q13 = q1(3);
q14 = q1(4);
q15 = q1(5);


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

T1_zeta1 = exp_se3(zeta_r*q11);
T1_zeta2 = exp_se3(zeta_r*q12);
T1_zeta3 = exp_se3(zeta_r*q13);
T1_zeta4 = exp_se3(zeta_p*q14);
T1_zeta5 = exp_se3(zeta_r*q15);

Te_seq(:,:,1) = T1_1*T1_zeta1*T1_2*T1_zeta2*T1_3*T1_zeta3*T1_4*T1_zeta4*T1_5*T1_zeta5*T1_p;


% UPS支链
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

    q = joint_q(:, i_limb);
    T_zeta1 = exp_se3(zeta_r*q(1));
    T_zeta2 = exp_se3(zeta_r*q(2));
    T_zeta3 = exp_se3(zeta_p*q(3));
    T_zeta4 = exp_se3(zeta_r*q(4));
    T_zeta5 = exp_se3(zeta_r*q(5));
    T_zeta6 = exp_se3(zeta_r*q(6));

    Te_seq(:,:,i_limb) = T01*T_zeta1*T12*T_zeta2*T23*T_zeta3*T34*T_zeta4*T45*T_zeta5*T56*T_zeta6*T_p;
end



end