function xi_seq = rebuild_xi_seq_from_p(p_seq)
% 根据当前 p_seq 重新计算各关节零位下的全局旋量坐标 xi_seq

    xi_seq = zeros(size(p_seq));

    zeta_r = [0;0;1;0;0;0];
    zeta_p = [0;0;0;0;0;1];

    % SPR 支链（1-6）
    T1 = zeros(4,4,6);
    for i_joint = 1 : 6
        T1(:,:,i_joint) = exp_se3(p_seq(:, i_joint));
    end
    xi_seq(:, 1) = adjoint_m(T1(:,:,1)) * zeta_r;
    xi_seq(:, 2) = adjoint_m(T1(:,:,1) * T1(:,:,2)) * zeta_r;
    xi_seq(:, 3) = adjoint_m(T1(:,:,1) * T1(:,:,2) * T1(:,:,3)) * zeta_r;
    xi_seq(:, 4) = adjoint_m(T1(:,:,1) * T1(:,:,2) * T1(:,:,3) * T1(:,:,4)) * zeta_p;
    xi_seq(:, 5) = adjoint_m(T1(:,:,1) * T1(:,:,2) * T1(:,:,3) * T1(:,:,4) * T1(:,:,5)) * zeta_r;

    % UPS 支链（7-34）
    for i_limb = 2 : 5
        T2 = zeros(4,4,7);
        for i_joint = 1 : 7
            T2(:,:,i_joint) = exp_se3(p_seq(:, (7*(i_limb-1) - 1) + i_joint));
        end

        xi_seq(:, 7*(i_limb-1) + 0) = adjoint_m(T2(:,:,1)) * zeta_r;
        xi_seq(:, 7*(i_limb-1) + 1) = adjoint_m(T2(:,:,1) * T2(:,:,2)) * zeta_r;
        xi_seq(:, 7*(i_limb-1) + 2) = adjoint_m(T2(:,:,1) * T2(:,:,2) * T2(:,:,3)) * zeta_p;
        xi_seq(:, 7*(i_limb-1) + 3) = adjoint_m(T2(:,:,1) * T2(:,:,2) * T2(:,:,3) * T2(:,:,4)) * zeta_r;
        xi_seq(:, 7*(i_limb-1) + 4) = adjoint_m(T2(:,:,1) * T2(:,:,2) * T2(:,:,3) * T2(:,:,4) * T2(:,:,5)) * zeta_r;
        xi_seq(:, 7*(i_limb-1) + 5) = adjoint_m(T2(:,:,1) * T2(:,:,2) * T2(:,:,3) * T2(:,:,4) * T2(:,:,5) * T2(:,:,6)) * zeta_r;
    end
end
