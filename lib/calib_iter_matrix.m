function [Jp_bar, Jp, Xi_bar] = calib_iter_matrix(joint_q, p_seq)
% 求对角块矩阵的过程可以用cell+blkdiag(cell{:})的方式写成循环

J = jacobian_space(joint_q, p_seq);  % 6*6 - 5
Omega = [zeros(3,3) eye(3); eye(3) zeros(3,3)];

Jp_blocks = cell(1, 5);
Xi_bar = zeros(6,6);

% 局部指数基公式
zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
zeta_p = [0;0;0;0;0;1];  % 平移基底

% SPR
T1 = zeros(4,4,6);
T1xi = zeros(4,4,5);
T1zeta = zeros(4,4,5);
U1_blocks = cell(1, 6);
U1_blocks{6} = eye(6);
D1_blocks = cell(1, 6);
D1_blocks{6} = eye(6);

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
        U1_blocks{i_joint} = adjoint_m(T_temp) * zeta_r;
    else
        T_zeta = exp_se3(zeta_p * joint_q(i_joint,1));
        U1_blocks{i_joint} = adjoint_m(T_temp) * zeta_p;
    end

    T1zeta(:,:,i_joint) = T_zeta;
    T1xi(:,:,i_joint) = T_temp * T_zeta / T_temp;
    D1_blocks{i_joint} = adjoint_m(T1xi(:,:,i_joint));
    
end

J1_passive = [J(:, 1:3, 1) J(:, 5, 1)];
Xi1 = null(J1_passive'*Omega);  % 对应论文公式2-42
Ap1 = blkdiag(Ap(p_seq(:,1)), ...
                Ap(p_seq(:,2)), ...
                Ap(p_seq(:,3)), ...
                Ap(p_seq(:,4)), ...
                Ap(p_seq(:,5)), ...
                Ap(p_seq(:,6)));

Lambda1 = blkdiag(eye(6), ...
                adj_m(T1(:,:,1)), ...
                adj_m(T1(:,:,1)*T1(:,:,2)), ...
                adj_m(T1(:,:,1)*T1(:,:,2)*T1(:,:,3)), ...
                adj_m(T1(:,:,1)*T1(:,:,2)*T1(:,:,3)*T1(:,:,4)), ...
                adj_m(T1(:,:,1)*T1(:,:,2)*T1(:,:,3)*T1(:,:,4)*T1(:,:,5))  );
Gamma1 = [eye(6), ...
        adj_m(T1xi(:,:,1)), ...
        adj_m(T1xi(:,:,1))*adj_m(T1xi(:,:,2)), ...
        adj_m(T1xi(:,:,1))*adj_m(T1xi(:,:,2))*adj_m(T1xi(:,:,3)), ...
        adj_m(T1xi(:,:,1))*adj_m(T1xi(:,:,2))*adj_m(T1xi(:,:,3))*adj_m(T1xi(:,:,4)), ...
        adj_m(T1xi(:,:,1))*adj_m(T1xi(:,:,2))*adj_m(T1xi(:,:,3))*adj_m(T1xi(:,:,4))*adj_m(T1xi(:,:,5))];

U1 = blkdiag(U1_blocks{:});
D1 = blkdiag(D1_blocks{:});

% 可验证与直接算Jp1相等，分开算是为了后续便于去除冗余参数
% Jp1 = [Ap(p_seq(:,1)), ...
%         adj_m(T1(:,:,1)*T1zeta(:,:,1)) * Ap(p_seq(:,2)), ...
%         adj_m(T1(:,:,1)*T1zeta(:,:,1)*T1(:,:,2)*T1zeta(:,:,2)) * Ap(p_seq(:,3)), ...
%         adj_m(T1(:,:,1)*T1zeta(:,:,1)*T1(:,:,2)*T1zeta(:,:,2)*T1(:,:,3)*T1zeta(:,:,3)) * Ap(p_seq(:,4)), ...
%         adj_m(T1(:,:,1)*T1zeta(:,:,1)*T1(:,:,2)*T1zeta(:,:,2)*T1(:,:,3)*T1zeta(:,:,3)*T1(:,:,4)*T1zeta(:,:,4)) * Ap(p_seq(:,5)), ...
%         adj_m(T1(:,:,1)*T1zeta(:,:,1)*T1(:,:,2)*T1zeta(:,:,2)*T1(:,:,3)*T1zeta(:,:,3)*T1(:,:,4)*T1zeta(:,:,4)*T1(:,:,5)*T1zeta(:,:,5)) * Ap(p_seq(:,6)) ];

Jp1 = Gamma1*Lambda1*Ap1;
Jp_blocks{1} = Xi1'*Omega*Jp1;  % 对应式2-45
Xi_bar(1:2, :) = Xi1'*Omega;


% UPS
T2 = zeros(4,4,7);
T2zeta = zeros(4,4,6);
T2xi = zeros(4,4,6);
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
        T2xi(:,:,i_joint) = T_temp * T_zeta / T_temp;
    end

    % 调试用
    J2_passive = [J(:, 1:2, i_limb) J(:, 4:6, i_limb)];
    % fprintf("i_limb = %d, rank(J2_passive) = %d\n",i_limb, rank(J2_passive));
    Xi2 = null(J2_passive'*Omega);

    Ap2 = blkdiag( ...
                Ap(p_seq(:, 7*(i_limb-1))), ...
                Ap(p_seq(:, 7*(i_limb-1)+1)), ...
                Ap(p_seq(:, 7*(i_limb-1)+2)), ...
                Ap(p_seq(:, 7*(i_limb-1)+3)), ...
                Ap(p_seq(:, 7*(i_limb-1)+4)), ...
                Ap(p_seq(:, 7*(i_limb-1)+5)), ...
                Ap(p_seq(:, 7*(i_limb-1)+6)) ...
                );

    Lambda2 = blkdiag( ...
                eye(6), ...
                adj_m(T2(:,:,1)), ...
                adj_m(T2(:,:,1)*T2(:,:,2)), ...
                adj_m(T2(:,:,1)*T2(:,:,2)*T2(:,:,3)), ...
                adj_m(T2(:,:,1)*T2(:,:,2)*T2(:,:,3)*T2(:,:,4)), ...
                adj_m(T2(:,:,1)*T2(:,:,2)*T2(:,:,3)*T2(:,:,4)*T2(:,:,5)), ...
                adj_m(T2(:,:,1)*T2(:,:,2)*T2(:,:,3)*T2(:,:,4)*T2(:,:,5)*T2(:,:,6)) ...
                    );

    Gamma2 = [ ...
        eye(6), ...
        adj_m(T2xi(:,:,1)), ...
        adj_m(T2xi(:,:,1))*adj_m(T2xi(:,:,2)), ...
        adj_m(T2xi(:,:,1))*adj_m(T2xi(:,:,2))*adj_m(T2xi(:,:,3)), ...
        adj_m(T2xi(:,:,1))*adj_m(T2xi(:,:,2))*adj_m(T2xi(:,:,3))*adj_m(T2xi(:,:,4)), ...
        adj_m(T2xi(:,:,1))*adj_m(T2xi(:,:,2))*adj_m(T2xi(:,:,3))*adj_m(T2xi(:,:,4))*adj_m(T2xi(:,:,5)), ...
        adj_m(T2xi(:,:,1))*adj_m(T2xi(:,:,2))*adj_m(T2xi(:,:,3))*adj_m(T2xi(:,:,4))*adj_m(T2xi(:,:,5))*adj_m(T2xi(:,:,6)) ...
            ];

    
    Jp2 = Gamma2*Lambda2*Ap2;
    Jp_blocks{i_limb} = Xi2'*Omega*Jp2;
    Xi_bar(i_limb+1, :) = Xi2' * Omega;
end

Jp = blkdiag(Jp_blocks{:});
Jp_bar = Xi_bar \ Jp;

end


function AdT = adj_m(T)
    R = T(1:3, 1:3);
    t = T(1:3, 4);
    AdT = [R zeros(3,3); skew(t)*R R];
end

function A = Ap(screw)
    omega = screw(1:3);
    nu = screw(4:6);
    
    Z = [skew(omega) zeros(3,3); skew(nu) skew(omega)];
    phi = norm(omega);  % 是旋转模长，不包括线速度

    if phi < 1e-5
        A = eye(6) + 0.5 * Z + (1/6) * (Z*Z);
    else
        S = sin(phi);
        C = cos(phi);
        
        c1 = (4 - phi*S - 4*C) / (2*phi^2);
        c2 = (4*phi - 5*S + phi*C) / (2*phi^3);
        c3 = (2 - phi*S - 2*C) / (2*phi^4);
        c4 = (2*phi - 3*S + phi*C) / (2*phi^5);
        
        Z2 = Z*Z;
        Z3 = Z2*Z;
        Z4 = Z3*Z;
        
        A = eye(6) + c1*Z + c2*Z2 + c3*Z3 + c4*Z4;
    end
end

function S = skew(w)
    S = [0,    -w(3),  w(2);
         w(3),  0,    -w(1);
        -w(2),  w(1),  0];
end