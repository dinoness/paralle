function J = jacobian_body(joint_q, p_seq)
% 计算物体雅可比

J = zeros(6,6,5);
T = keni_sol_forward_once(joint_q, p_seq);

J_space = jacobian_space(joint_q, p_seq);

for i_limb = 1 : 5
    if i_limb ~= 1
        n_joint = 6;
    else
        n_joint = 5;
    end

    for i_joint = 1 : n_joint
        J(:, i_joint, i_limb) = adjoint(trans_inv(T(:,:,i_limb)),J_space(:, i_joint, i_limb));
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

function T_inv = trans_inv(T)
    R = T(1:3, 1:3);
    t = T(1:3, 4);
    T_inv = [R' -R'*t; 0 0 0 1];    
end