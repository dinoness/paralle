function R = exp_so3(phi)
%  EXP_SO3 so(3)李代数指数映射
%  input: xi 3*1向量，旋转
%  output: T 3*3矩阵
theta = norm(phi);
if theta < 1e-10
        R = eye(3);
        return;
end
phi_unit = phi / theta;
phi_up = [           0  -phi_unit(3)  phi_unit(2);
           phi_unit(3)             0 -phi_unit(1);
          -phi_unit(2)  phi_unit(1)             0]; % 反对称矩阵
R = cos(theta)*eye(3) + (1-cos(theta))* (phi_unit * phi_unit') + sin(theta)*phi_up;

end