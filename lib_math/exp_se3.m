function T = exp_se3(xi)
%  EXP_SE3 se(3)李代数指数映射
%  input: xi 6*1向量，旋转/位移
%  output: T 4*4矩阵
phi = xi(1:3);
rou = xi(4:6);
o13 = zeros(1,3);

theta = norm(phi);
if theta < 1e-10
    R = eye(3);
    T = [R rou; o13 1];
    return;
end
phi_unit = phi / theta;
phi_up = [           0  -phi_unit(3)  phi_unit(2);
           phi_unit(3)             0 -phi_unit(1);
          -phi_unit(2)  phi_unit(1)             0]; % 反对称矩阵

R = cos(theta)*eye(3) + (1-cos(theta))* (phi_unit * phi_unit') + sin(theta)*phi_up;

J = sin(theta) / theta * eye(3) ...
    + (1 - sin(theta) / theta) * (phi_unit * phi_unit') ...
    + (1 - cos(theta)) / theta * phi_up;
T = [R J*rou; o13 1];

% fprintf("theta = %f\n", theta)
% fprintf("phi = \n")
% disp(phi)
% fprintf("J = \n")
% disp(J)

end


% ===========test data============
% phi = pi/4/sqrt(6)*[1;-2;-1];
% rou = [0;0;0.5];
% xi = [phi;rou];

% ===========paper data============
% p1 = [0;0;2.0944;0.2056;0.3560;0];
% p2 = [-1.5708;0;0;0;0;0];
% p3 = [0;0;0;0;0.5;0];
% p4 = [1.5708;0;0;-0.15;0;0];
% p5 = [0;1.5708;0;0;0;0];
% p6 = [1.5835;-0.9142;-1.5835;-0.0325;-0.1343;-0.1];

% function T = exp_se3(xi)
% % se3_exp    Compute the exponential map from se(3) to SE(3).
% %   T = se3_exp(xi) calculates the matrix exponential of a twist xi in se(3).
% %
% %   Input:
% %       xi - 6-element vector [omega; v]
% %            omega : 3x1 rotation vector (axis * angle)
% %            v     : 3x1 translation vector
% %   Output:
% %       T  - 4x4 homogeneous transformation matrix in SE(3):
% %            [R, t; 0 0 0 1]
% 
%     % Extract rotation and translation parts
%     omega = xi(1:3);
%     v     = xi(4:6);
% 
%     theta = norm(omega);            % rotation angle
%     tol   = 1e-10;                  % threshold for small angle
% 
%     % Pre‑compute skew‑symmetric matrix and its square
%     skew_w  = skew(omega);
%     skew_w2 = skew_w * skew_w;
% 
%     % Compute coefficients for the exponential formula
%     if theta < tol
%         % Use Taylor series expansions for numerical stability
%         % sin(theta)/theta = 1 - theta^2/6 + theta^4/120 - ...
%         s_over_theta = 1 - theta^2/6 + theta^4/120;
%         % (1 - cos(theta))/theta^2 = 1/2 - theta^2/24 + theta^4/720 - ...
%         one_minus_cos_over_theta2 = 0.5 - theta^2/24 + theta^4/720;
%         % (theta - sin(theta))/theta^3 = 1/6 - theta^2/120 + theta^4/5040 - ...
%         theta_minus_sin_over_theta3 = 1/6 - theta^2/120 + theta^4/5040;
%     else
%         % Direct evaluation (well‑conditioned)
%         s_over_theta = sin(theta) / theta;
%         one_minus_cos_over_theta2 = (1 - cos(theta)) / theta^2;
%         theta_minus_sin_over_theta3 = (theta - sin(theta)) / theta^3;
%     end
% 
%     % Rodrigues formula for rotation matrix
%     R = eye(3) + s_over_theta * skew_w + one_minus_cos_over_theta2 * skew_w2;
% 
%     % Translation part: t = V * v, where V is given by
%     V = eye(3) + one_minus_cos_over_theta2 * skew_w + theta_minus_sin_over_theta3 * skew_w2;
%     t = V * v;
% 
%     % Assemble homogeneous transformation matrix
%     T = [R, t; 0 0 0 1];
% end
% 
% function S = skew(w)
% % skew    Returns the 3x3 skew‑symmetric matrix of a 3‑vector.
% %   S = skew(w) generates the matrix such that S * v = w × v for any v.
%     S = [0,    -w(3),  w(2);
%          w(3),  0,    -w(1);
%         -w(2),  w(1),  0];
% end