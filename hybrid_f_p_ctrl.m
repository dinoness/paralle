%% 待解决问题
% 1. 不同运动模式的切换
% 2. 力分量的控制律
% 3. 运动学正解的可行性
% 4. 力传感器的选型、信号接入与滤波处理


%% 输入
% 期望力+期望姿态
f_seq = [10, 10];
p_seq = [0, 1;
         0, -1;
         -800, -801;
         0, 5;
         2, 3];

% 初始状态
f0 = 9;
p0 = [0;0;-799;0;0];

p_last = p0;
f_last = f0;
%% 控制周期
n_seq = length(f_seq);

for i = 1 : n_seq
    p_ref = p_seq(:, i);
    f_ref = f_seq(i);

    p_cur = p_last; % 正解获得，否则用上一个时刻的位姿
    f_cur = f_last; % 读取传感器数据获得，否则用上一个时刻的期望力

    % =======力控制器=======
    % 力误差
    f_err = f_ref - f_cur;

    % 力的控制律，输出位移，待完善
    dz = 0.01 * f_err;

    % z轴分量
    % 刀具姿态有3种选择，1正解获得的当前姿态，2上一个周期的期望姿态，3本周期的期望姿态
    phi = p_cur(4);
    theta = p_cur(5);
    R_z = [cos(phi)*sin(theta); sin(phi)*sin(theta); cos(theta)];

    % 将位移进行变换
    p_force_ctrl = [dz*R_z; 0; 0];


    % =======位置控制器=======
    % 相对位移量，此处用矩阵，实际操作时直接取即可
    S = diag([1 1 0 1 1]);

    % 直接做差，暂不进行特殊控制，结束
    p_position_ctrl = S * (p_ref - p_cur);

    % =======输入控制器的位移量=======
    p_out = p_position_ctrl + p_force_ctrl;

    p_last = p_cur;
    f_last = f_cur;
end