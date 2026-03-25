function [T0_seq, Txi_seq] = transform_matrix_cal(p_seq, joint_q)
% 从结构参数和关节量计算转移矩阵
% 待写=====================

T0_seq = cell(1, 5);
Txi_seq = cell(1, 5);


% 局部指数基公式
zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
zeta_p = [0;0;0;0;0;1];  % 平移基底



end