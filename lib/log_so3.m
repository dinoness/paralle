function [phi, theta] = log_so3(R)
%  LOG_SO3 so(3)李群对数映射
%  input: R 3*3矩阵
%  output: phi 3*1向量，旋转
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

end