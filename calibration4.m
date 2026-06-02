% 参考文献：c02 并联机构的运动学误差建模及参数可辨识性分析_孔令雨
% L-M求解
clear
path_add();
fprintf('>>>= start (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));


%% 参数集
% -----struct-parameter------
basic_paras = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'm');  % 单位意思是程序中用到的单位
l_max = basic_paras.l_max;
l_min = basic_paras.l_min;
R1 = basic_paras.R1;
R2 = basic_paras.R2;
H = basic_paras.H;
r1 = basic_paras.r1;
r2 = basic_paras.r2;
h = basic_paras.h;
L_tool = basic_paras.L_tool;
limb_dir = basic_paras.limb_dir;
B = basic_paras.B;
P_m = basic_paras.P_m;
l0_seq = basic_paras.l0_seq;
joint_u_angle_tilt = basic_paras.joint_u_angle_tilt;
unit_para = basic_paras.unit_para;
% -----end-struct-parameter------

err_max = 3e-6;
keni_forward_err_max = 1e-8;
loop_max = 100;
para_err_std = 0.0010;  % 制造误差标准差（m/rad），模拟真实参数与名义参数的偏差
a_dis = 1e-5;

% ----- input data ------
% 参考序列生成
x_seq = [-100 0 100] * unit_para;
y_seq = [-100 0 100] * unit_para;
z_seq = [-900 -1000] * unit_para;
phi_seq = deg2rad(-120 : 120 : 120); 
theta_seq = deg2rad([0 5 10]);
Pos_ref_seq = zeros(5, 3*3*3*3*2);
for ix = 1:length(x_seq)
    for iy = 1:length(y_seq)
        for iz = 1:length(z_seq)
            for iphi = 1:length(phi_seq)
                for itheta = 1:length(theta_seq)
                    Pos_ref_seq(:, (ix-1)*54+(iy-1)*18+(iz-1)*9+(iphi-1)*3+(itheta)) = [x_seq(ix); y_seq(iy); z_seq(iz); phi_seq(iphi); theta_seq(itheta)];
                end
            end
        end
    end
end

% x_seq = [-100 0 100] * unit_para;
% y_seq = [-100 0 100] * unit_para;
% z_seq = [-800 -850 -900] * unit_para;
% theta_seq = deg2rad([0 5 10]);
% Pos_ref_seq = zeros(5, 3*3*3*3);
% for ix = 1:3
%     for iy = 1:3
%         for iz = 1:3
%             for itheta = 1:3
%                 Pos_ref_seq(:, (ix-1)*27+(iy-1)*9+(iz-1)*3+(itheta)) = [x_seq(ix); y_seq(iy); z_seq(iz); 0; theta_seq(itheta)];
%             end
%         end
%     end
% end
% 仿真测量数据通过以下流程生成：
%   名义逆解(关节指令) → 真实模型正解(含制造误差) → 输出位姿(模拟测量值)
% 与直接对位姿加扰动不同，此方法更真实地反映参数误差在位姿空间中的传播
% Pos_m_seq = [0.01;-0.01;-600.02;0.001;-0.001];  % line=5 colum=n
% Pos_ref_seq = [0;30;-600;0;10];  % line=5 colum=n  角度的单位是° **一列为一组**
seq_len = length(Pos_ref_seq(1, :));
Pos_err_seq = zeros(5, seq_len);  % 位姿估计误差，优化的目标
Pos_delta_seq = zeros(5, seq_len);  % 位姿扰动序列
% ----- end input data ------

%% 标定步骤
[p_seq_nom, xi_seq_nom] = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

% 创建"真实"模型：在名义参数上叠加制造误差
rng(0313);
p_seq_true = p_seq_nom + para_err_std * randn(6, 34);

T_cal_seq = zeros(4, 4, seq_len);
T_measure_seq = zeros(4, 4, seq_len);
joint_seq_iter = zeros(6, 5, seq_len);
err_seq_iter = zeros(6*seq_len, 1);
p_seq_iter = p_seq_nom;  % 结构参数序列
xi_seq_iter = xi_seq_nom;

%% 仿真测量数据生成
%   步骤1：名义逆解 → 主动关节量（"指令值"）
%   步骤2：真实模型正解 → 实际输出位姿（"测量值"）
for im = 1 : seq_len
    T_cal_seq(:,:,im) = pos2trans(Pos_ref_seq(:, im), B,'unit','rad');

    % 步骤1：名义模型逆解得到所有关节量（含主动+被动）
    joint_q_ref = keni_sol_inverse(T_cal_seq(:,:,im), B, l0_seq, P_m, p_seq_nom);

    % 步骤2：固定主动关节（移动副），用真实模型正解得到实际位姿
    [T_measure_seq(:,:,im), ~] = keni_sol_forward(joint_q_ref, p_seq_true, keni_forward_err_max);
    % rng(0313+im);  % 随机数种子
    % screw_temp = log_se3(T_cal_seq(:,:,im)) + a_dis * rand(6, 1);
    % T_measure_seq(:,:,im) = exp_se3(screw_temp);  % 通过添加扰动获得实际位姿（之后用数据替代

    % 标定初始关节量（名义逆解值）
    joint_seq_iter(:,:,im) = joint_q_ref;

    % 初始位姿误差
    err_seq_iter(6*(im-1)+1 : 6*im) = log_se3(T_measure_seq(:,:,im) / T_cal_seq(:,:,im));

end


%% 引入真实位姿，判断误差是否小于容差，是则结束，否则进入迭代
calib_loop = 0;
Jp_bar = zeros(6*seq_len, 112);
err_list = zeros(loop_max+1,1);
lambda = 1e-1;

err_cur = norm(err_seq_iter);
err_list(1) = err_cur;
while err_cur > err_max
    [U, N, V_prep, M] = calib_iter_row_space_matrix(xi_seq_iter);
    % size(U{1})
    % size(V_prep{1})
    % size(V{1})
    [Lambda, Ap] = calib_iter_restore_matrix(p_seq_iter);
    % 更新运动学参数
    for im = 1 : seq_len
        % Jp_bar(6*(im-1)+1 : 6*(im-1) + 6, :) = calib_iter_matrix(joint_seq_iter(:,:,im), p_seq_iter);
        Jp_bar(6*(im-1)+1 : 6*(im-1) + 6, :) = calib_iter_matrix2(joint_seq_iter(:,:,im), p_seq_iter, xi_seq_iter, U, V_prep);
    end
    

    pk = (Jp_bar' * Jp_bar + lambda * eye(size(Jp_bar, 2))) \ (Jp_bar' * err_seq_iter);
    
    delta_p_seq = restore_full_param_increment(pk, U, V_prep, Lambda, Ap);  % 依据式(3-28)恢复为原始运动学参数
    p_seq_vec = p_seq_iter(:) + delta_p_seq(:);  % A(:)矩阵按列排列成列向量
    p_seq_iter = reshape(p_seq_vec, size(p_seq_iter, 1), size(p_seq_iter, 2));  % 将向量重排为矩阵
    % p_seq_vec = p_seq_iter(:) + 0.001 * pinv(Jp_bar' * Jp_bar) * Jp_bar' * err_seq_iter;
    

    xi_seq_iter = rebuild_xi_seq_from_p(p_seq_iter);  % 更新零位全局旋量坐标，供下一轮迭代使用
    
    for im = 1 : seq_len
        % 新运动学参数下的正解
        [T_cal_seq(:,:,im), joint_seq_iter(:,:,im)] = keni_sol_forward(joint_seq_iter(:,:,im), p_seq_iter, keni_forward_err_max);
        % 计算位姿误差
        err_seq_iter(6*(im-1)+1 : 6*im) = log_se3(T_measure_seq(:,:,im) / T_cal_seq(:,:,im));
    end

    
    % 阻尼系数迭代
    err_new = norm(err_seq_iter);
    delta_err = err_new - err_cur;  % 实际下降
    delta_qk = -1 *( (Jp_bar'*err_seq_iter)'*pk + 0.5*pk'*(Jp_bar'*Jp_bar)*pk );  % 期望下降
    eta = delta_err / delta_qk;
            
    % 系数可调整
    if eta > 0.75
        lambda = lambda * 0.5;
    elseif eta < 0.25
        lambda = lambda * 3;
    end
    

    calib_loop = calib_loop + 1;
    err_cur = err_new;
    err_list(calib_loop+1) = err_cur;
    if calib_loop > loop_max
        break;
    end
    
    if(rem(calib_loop, 10) == 0)
        fprintf("loop = %d\n", calib_loop);
    end
end

fig = figure('Color', [1 1 1]);
plot(0:calib_loop, err_list(1:calib_loop+1),'linewidth',1.5)
set(gca, 'YScale', 'log');  % 对数坐标轴
grid on
set(gca,'linewidth',1.5,'fontsize',15,'fontname','Times New Roman');
set(gcf,'unit','centimeters','position',[10 10 14 8]);
xlabel('迭代次数', 'FontSize', 14, 'FontName', '微软雅黑', 'FontWeight', 'bold', 'Color', 'black');
ylabel('残余误差', 'FontSize', 14, 'FontName', '微软雅黑', 'FontWeight', 'bold', 'Color', 'black');


% keni_sol_forward_once(joint_seq_iter(:,:,1),p_seq_iter)
% T_measure_seq(:,:,1)
% T_cal_seq(:,:,1)
