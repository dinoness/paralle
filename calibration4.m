% 参考文献：c02 并联机构的运动学误差建模及参数可辨识性分析_孔令雨
clear
addpath(genpath('./lib'));
%% 参数集
%--------struct parameter--------
T = readtable('parameters.xlsx', 'Range', 'A2:B12');
paras = table2array(T(:, 2));
l_max = paras(1);
l_min = paras(2);  % 670
R1 = paras(3);  % 550
R2 = paras(4);  % 500
H = paras(5);  % 0
r1 = paras(6);  % 100
r2 = paras(7);  % 80
h = paras(8);  % 10

limb_dir = [pi/2; 7*pi/6; -pi/6; pi/6; 5*pi/6];
B1 = [R1*cos(pi/2);   R1*sin(pi/2);   0];
B2 = [R1*cos(7*pi/6); R1*sin(7*pi/6); 0];
B3 = [R1*cos(-pi/6);  R1*sin(-pi/6);  0];
B4 = [R2*cos(pi/6);   R2*sin(pi/6);   H];
B5 = [R2*cos(5*pi/6); R2*sin(5*pi/6); H];
B = [B1 B2 B3 B4 B5];

% move plant parameter
P1_m = [r1*cos(pi/2);   r1*sin(pi/2);   0];
P2_m = [r1*cos(7*pi/6); r1*sin(7*pi/6); 0];
P3_m = [r1*cos(-pi/6);  r1*sin(-pi/6);  0];
P4_m = [r2*cos(pi/6);   r2*sin(pi/6);   h];
P5_m = [r2*cos(5*pi/6); r2*sin(5*pi/6); h];
P_m = [P1_m P2_m P3_m P4_m P5_m];
P_v = zeros(3, 5);  % 只变换了方向，没变换起点
P = zeros(3, 5);    % 末端点坐标


l0 = 600;  % 默认初始支链长度
l0_seq = [l0;l0;l0;l0;l0];
joint_u_angle_tilt = 155 / 180 * pi;
% -----end-struct-parameter------

err_max = 1e-3;
loop_max = 5;
a_dis = 0.1;  % 扰动幅度

% ----- input data ------
% 参考序列生成
x_seq = [-5 0 5];
y_seq = [-5 0 5];
z_seq = [-600 -615 -630];
theta_seq = [-5 0 5];
Pos_ref_seq = zeros(5, 3*3*3);
for ix = 1:3
    for iy = 1:3
        for iz = 1:3
            % for itheta = 1:3
            %     Pos_ref_seq(:, (ix-1)*27+(iy-1)*9+(iz-1)*3+(itheta)) = [x_seq(ix); y_seq(iy); z_seq(iz); 0; theta_seq(itheta)];
            % end
            Pos_ref_seq(:, (ix-1)*9+(iy-1)*3+iz) = [x_seq(ix); y_seq(iy); z_seq(iz); 0; 0];
        end
    end
end
% 当前难题，不知道测量得到的数据格式是什么，所以这里假定能通过数据处理软件得到三维位姿，此处用旋量扰动代替
% Pos_m_seq = [0.01;-0.01;-600.02;0.001;-0.001];  % line=5 colum=n
% Pos_ref_seq = [0;30;-600;0;10];  % line=5 colum=n  角度的单位是° **一列为一组**
seq_len = length(Pos_ref_seq(1, :));
Pos_err_seq = zeros(5, seq_len);  % 位姿估计误差，优化的目标
Pos_delta_seq = zeros(5, seq_len);  % 位姿扰动序列
% ----- end input data ------

%% 标定步骤
p_seq_nom = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);
T_cal_seq = zeros(4, 4, seq_len);
T_measure_seq = zeros(4, 4, seq_len);
joint_seq_iter = zeros(6, 5, seq_len);
err_seq_iter = zeros(6*seq_len, 1);
p_seq_iter = p_seq_nom;  % 结构参数序列


%% 计算名义参数下，理想位姿的运动学正逆解
for im = 1 : seq_len
    % 名义与测量位姿
    T_cal_seq(:,:,im) = pos2trans(Pos_ref_seq(:, im), B);
    screw_temp = log_se3(T_cal_seq(:,:,im)) + a_dis * rand(6, 1);
    T_measure_seq(:,:,im) = exp_se3(screw_temp);  % 通过添加扰动获得实际位姿（之后用数据替代）

    % 初始关节量
    joint_q_ref = keni_sol_inverse(T_cal_seq(:,:,im), B, l0_seq, P_m, p_seq_nom);
    joint_seq_iter(:,:,im) = joint_q_ref;

    % 计算测量值与参考值之间的误差
    err_seq_iter(6*(im-1)+1 : 6*im) = log_se3(T_measure_seq(:,:,im) / T_cal_seq(:,:,im));
end


%% 引入真实位姿，判断误差是否小于容差，是则结束，否则进入迭代
calib_loop = 0;
Jp_bar = zeros(6*seq_len, 204);
err_list = zeros(loop_max,1);
while norm(err_seq_iter) > err_max
    % 更新运动学参数
    for im = 1 : seq_len
        Jp_bar(6*(im-1)+1 : 6*(im-1) + 6, :) = calib_iter_matrix(joint_seq_iter(:,:,im), p_seq_iter);
    end

    
    p_seq_vec = p_seq_iter(:) + 0.000001 * pinv(Jp_bar' * Jp_bar) * Jp_bar' * err_seq_iter;  % A(:)矩阵按列排列成列向量
    p_seq_iter = reshape(p_seq_vec, size(p_seq_iter, 1), size(p_seq_iter, 2));  % 将向量重排为矩阵
    
    for im = 1 : seq_len
        % 新运动学参数下的正解
        T_cal_seq(:,:,im) = keni_sol_forward(joint_seq_iter(:,:,im), p_seq_iter, err_max);
        % 计算位姿误差
        err_seq_iter(6*(im-1)+1 : 6*im) = log_se3(T_measure_seq(:,:,im) / T_cal_seq(:,:,im));
    end

    calib_loop = calib_loop + 1;
    err_list(calib_loop) = norm(err_seq_iter);
    if calib_loop > loop_max
        break;
    end
end

fig = figure('Color', [1 1 1]);
plot(err_list(1:calib_loop))


% 当前的正解容易发散，当结构参数变化过大，如其中一个量变了10
% 位姿迭代矩阵步幅大，且迭代方向不对，即使补偿乘了1e-6的系数，误差也一直增大