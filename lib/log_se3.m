function xi = log_se3(T)
%  EXP_SE3 se(3)李群对数映射
%  input: T 4*4矩阵
%  output: xi 6*1向量，旋转/位移
R = T(1:3, 1:3);
t = T(1:3, 4);

% theta = acos((trace(R) - 1)/2);
% [V, D] = eig(R);  %
% if abs(real(D(3,3)) - 1) < 0.00001
%     phi_unit = V(:, 3);
%     phi = theta * V(:, 3);
% else
%     fprintf("R:\n")
%     disp(R);
%     fprintf("V:\n")
%     disp(V);
%     fprintf("D:\n")
%     disp(D);
%     returnw
% end
theta = acos((trace(R) - 1)/2);
axis = 1 / (2 * sin(theta)) * [R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)];
phi = theta * axis;

if theta < 1e-10
    return;
end

% phi_up = [           0  -phi_unit(3)  phi_unit(2);
%            phi_unit(3)             0 -phi_unit(1);
%           -phi_unit(2)  phi_unit(1)             0]; % 反对称矩阵
% 
% J = sin(theta) / theta * eye(3) ...
%     + (1 - sin(theta) / theta) * (phi_unit * phi_unit') ...
%     + (1 - cos(theta)) / theta * phi_up;

J = sin(theta) / theta * eye(3) ...
    + (1 - sin(theta) / theta) * (axis * axis') ...
    + (1 - cos(theta)) / theta * skew(axis);

rou = J \ t;
xi = [phi;rou];

% fprintf("theta = %f\n", theta)
% fprintf("phi = \n")
% disp(phi)
% fprintf("J = \n")
% disp(J)


end


function S = skew(w)
% skew    Returns the 3x3 skew‑symmetric matrix of a 3‑vector.
%   S = skew(w) generates the matrix such that S * v = w × v for any v.
    S = [0,    -w(3),  w(2);
         w(3),  0,    -w(1);
        -w(2),  w(1),  0];
end