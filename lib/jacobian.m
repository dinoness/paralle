function J = jacobian(joint_q, p_seq)
% 计算速度雅克比

J = zeros(6,6,5);

% 局部指数基公式
zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
zeta_p = [0;0;0;0;0;1];  % 平移基底

% SPR
T1 = zeros(4,4,6);
T1xi = zeros(4,4,5);
T1zeta = zeros(4,4,5);

for i_joint = 1 : 6
    T1(:,:,i_joint) = exp_se3(p_seq(:, i_joint));
end

% 计算转移矩阵
for i_joint = 1 : 5
    T_temp = eye(4);
    for i_T = 1 : i_joint
        T_temp = T_temp * T1(:,:,i_T);
    end

    if i_joint ~= 4
        T_zeta = exp_se3(zeta_r * joint_q(i_joint,1));
    else
        T_zeta = exp_se3(zeta_p * joint_q(i_joint,1));
    end

    T1zeta(:,:,i_joint) = T_zeta;
    T1xi(:,:,i_joint) = T_temp * T_zeta / T_temp;
end

% 计算雅可比
for i_joint = 1 : 5

    if i_joint > 1
        T_temp2 = T_temp2 * T1zeta(:,:,i_joint-1) * T1(:,:,i_joint);
    else
        T_temp2 = T1(:,:,i_joint);
    end

    if i_joint ~= 4
        J(:,i_joint,1) = adjoint(T_temp2, zeta_r);
    else
        J(:,i_joint,1) = adjoint(T_temp2, zeta_p);
    end
end


% 全局法，旋量本身要零位的，而算伴随矩阵要带上现有位姿的
% xi1 = adjoint(T1(:,:,1),zeta_r);
% xi2 = adjoint(T1(:,:,1)*T1(:,:,2),zeta_r);
% xi3 = adjoint(T1(:,:,1)*T1(:,:,2)*T1(:,:,3),zeta_r);
% xi4 = adjoint(T1(:,:,1)*T1(:,:,2)*T1(:,:,3)*T1(:,:,4),zeta_p);
% xi5 = adjoint(T1(:,:,1)*T1(:,:,2)*T1(:,:,3)*T1(:,:,4)*T1(:,:,5),zeta_r);

% J1 = xi1;
% J2 = adjoint(T1xi(:,:,1), xi2);
% J3 = adjoint(T1xi(:,:,1)*T1xi(:,:,2), xi3);
% J4 = adjoint(T1xi(:,:,1)*T1xi(:,:,2)*T1xi(:,:,3), xi4);
% J5 = adjoint(T1xi(:,:,1)*T1xi(:,:,2)*T1xi(:,:,3)*T1xi(:,:,4), xi5);
% [J1 J2 J3 J4 J5]

% UPS
T2 = zeros(4,4,7);
T2zeta = zeros(4,4,6);
for i_limb = 2 : 5
    for i_joint = 1 : 7
        T2(:,:,i_joint) = exp_se3(p_seq(:, (7*(i_limb-1) - 1) + i_joint));
    end

    for i_joint = 1 : 6
        T_temp = eye(4);
        for i_T = 1 : i_joint
            T_temp = T_temp * T2(:,:,i_T);
        end

        if i_joint ~= 3
            T_zeta = exp_se3(zeta_r * joint_q(i_joint,i_limb));
        else
            T_zeta = exp_se3(zeta_p * joint_q(i_joint,i_limb));
        end

        T2zeta(:,:,i_joint) = T_zeta;
    end

    for i_joint = 1 : 6

        if i_joint > 1
            T_temp2 = T_temp2 * T2zeta(:,:,i_joint-1) * T2(:,:,i_joint);
        else
            T_temp2 = T2(:,:,i_joint);
        end

        if i_joint ~= 3
            J(:,i_joint,i_limb) = adjoint(T_temp2, zeta_r);
        else
            J(:,i_joint,i_limb) = adjoint(T_temp2, zeta_p);
        end
    end
end


end


function j = adjoint(T, xi)
    R = T(1:3, 1:3);
    t = T(1:3, 4);
    j = [R zeros(3,3); skew(t)*R R] * xi;
end


function S = skew(w)
    S = [0,    -w(3),  w(2);
         w(3),  0,    -w(1);
        -w(2),  w(1),  0];
end