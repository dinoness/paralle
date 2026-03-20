function [out1, out2] = t_elements(joint_q, p_seq)
J = jacobian_space(joint_q, p_seq);  % 6*6 - 5
J2 = J(:,:,2);
J2_passive = [J2(:, 1:2) J2(:,4:6)];
Omega = [zeros(3,3) eye(3); eye(3) zeros(3,3)];
out1 = null(J2_passive'*Omega);

out2 = blkdiag(Ap(p_seq(:,7)), Ap(p_seq(:,8)), Ap(p_seq(:,9)), Ap(p_seq(:,10)), Ap(p_seq(:,11)), Ap(p_seq(:,12)), Ap(p_seq(:,13)));
end


function A = Ap(screw)
    omega = screw(1:3);
    nu = screw(4:6);
    
    Z = [skew(omega) zeros(3,3); skew(nu) skew(omega)];
    phi = norm(screw);

    if phi == 0
        A = eye(6) + 0.5 * Z;
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