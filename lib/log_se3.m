function xi = log_se3(T)
%  EXP_SE3 se(3)李群对数映射
%  input: T 4*4矩阵
%  output: xi 6*1向量，旋转/位移
R = T(1:3, 1:3);
t = T(1:3, 4);

theta = acos((trace(R) - 1)/2);
[V, D] = eig(R);  %
if abs(real(D(3,3)) - 1) < 0.00001
    phi = theta * V(:, 3);
else
    fprintf("R:\n")
    disp(R);
    fprintf("V:\n")
    disp(V);
    fprintf("D:\n")
    disp(D);
end

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

end