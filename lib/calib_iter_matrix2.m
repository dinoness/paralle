function [Jp_bar, Jp, Xi_bar] = calib_iter_matrix2(joint_q, p_seq, xi_seq, U, V_prep)
% 求对角块矩阵的过程可以用cell+blkdiag(cell{:})的方式写成循环

J = jacobian_space(joint_q, p_seq);  % 6*6 - 5
Omega = [zeros(3,3) eye(3); eye(3) zeros(3,3)];

Jp_blocks = cell(1, 5);
Xi_bar = zeros(6,6);


% 局部指数基公式
zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
zeta_p = [0;0;0;0;0;1];  % 平移基底

% SPR
T1xi = zeros(4,4,5);
D1_blocks = cell(1, 6);
D1_blocks{6} = eye(6);


T1 = zeros(4,4,6);
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

    T1xi(:,:,i_joint) = T_temp * T_zeta / T_temp;
    D1_blocks{i_joint} = eye(6) - adjoint_m(T1xi(:,:,i_joint));
    
end

J1_passive = [J(:, 1:3, 1) J(:, 5, 1)];
Xi1 = null(J1_passive'*Omega);  % 对应论文公式2-42

Gamma1 = [eye(6), ...
        adjoint_m(T1xi(:,:,1)), ...
        adjoint_m(T1xi(:,:,1))*adjoint_m(T1xi(:,:,2)), ...
        adjoint_m(T1xi(:,:,1))*adjoint_m(T1xi(:,:,2))*adjoint_m(T1xi(:,:,3)), ...
        adjoint_m(T1xi(:,:,1))*adjoint_m(T1xi(:,:,2))*adjoint_m(T1xi(:,:,3))*adjoint_m(T1xi(:,:,4)), ...
        adjoint_m(T1xi(:,:,1))*adjoint_m(T1xi(:,:,2))*adjoint_m(T1xi(:,:,3))*adjoint_m(T1xi(:,:,4))*adjoint_m(T1xi(:,:,5))];

D1 = blkdiag(D1_blocks{:});


Jp1 = Gamma1*D1*U{1}*V_prep{1};
Jp_blocks{1} = Xi1'*Omega*Jp1;  % 对应式2-45
Xi_bar(1:2, :) = Xi1'*Omega;


% UPS
T2 = zeros(4,4,7);
T2xi = zeros(4,4,6);
D2_blocks = cell(1, 7);
D2_blocks{7} = eye(6);
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

        T2xi(:,:,i_joint) = T_temp * T_zeta / T_temp;
        D2_blocks{i_joint} = eye(6) - adjoint_m(T2xi(:,:,i_joint));
    end

    % 调试用
    J2_passive = [J(:, 1:2, i_limb) J(:, 4:6, i_limb)];
    % fprintf("i_limb = %d, rank(J2_passive) = %d\n",i_limb, rank(J2_passive));
    Xi2 = null(J2_passive'*Omega);


    Gamma2 = [ ...
        eye(6), ...
        adjoint_m(T2xi(:,:,1)), ...
        adjoint_m(T2xi(:,:,1))*adjoint_m(T2xi(:,:,2)), ...
        adjoint_m(T2xi(:,:,1))*adjoint_m(T2xi(:,:,2))*adjoint_m(T2xi(:,:,3)), ...
        adjoint_m(T2xi(:,:,1))*adjoint_m(T2xi(:,:,2))*adjoint_m(T2xi(:,:,3))*adjoint_m(T2xi(:,:,4)), ...
        adjoint_m(T2xi(:,:,1))*adjoint_m(T2xi(:,:,2))*adjoint_m(T2xi(:,:,3))*adjoint_m(T2xi(:,:,4))*adjoint_m(T2xi(:,:,5)), ...
        adjoint_m(T2xi(:,:,1))*adjoint_m(T2xi(:,:,2))*adjoint_m(T2xi(:,:,3))*adjoint_m(T2xi(:,:,4))*adjoint_m(T2xi(:,:,5))*adjoint_m(T2xi(:,:,6)) ...
            ];
    
    D2 = blkdiag(D2_blocks{:});

    
    Jp2 = Gamma2*D2*U{i_limb}*V_prep{i_limb};
    Jp_blocks{i_limb} = Xi2'*Omega*Jp2;
    Xi_bar(i_limb+1, :) = Xi2' * Omega;
end

Jp = blkdiag(Jp_blocks{:});
Jp_bar = Xi_bar \ Jp;

end
