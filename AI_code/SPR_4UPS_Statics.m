% =========================================================================
% 主函数：SPR_4UPS_Statics
% 功能：计算给定外力和位姿下，SPR-4UPS机构的5个主动关节驱动力与1个约束力
% 输入：
%   T_mat: 动平台的4x4齐次变换矩阵 (当前位姿)
%   W_ext: 6x1 列向量，施加在动平台中心的外部力旋量 [Fx; Fy; Fz; Mx; My; Mz]
%   params: 结构体，包含机构的几何参数 (基座铰链点B，动平台局部铰链点A_local)
%   R_axis_global: (可选) SPR支链中R关节在全局坐标系下的方向向量，用于计算约束旋量
% 输出：
%   tau_active: 5x1 列向量，5个P关节的主动驱动力
%   f_constraint: 标量，SPR支链提供的约束力/力矩大小
%   G_matrix: 6x6 总体力雅可比矩阵
% =========================================================================
function [tau_active, f_constraint, G_matrix] = SPR_4UPS_Statics(T_mat, W_ext, params)
    
    % 解析动平台位姿
    R_mat = T_mat(1:3, 1:3); % 旋转矩阵
    P_pos = T_mat(1:3, 4);   % 位置向量

    % 提取铰链点参数
    B = params.B;                 % 3x5 矩阵，基座铰链点坐标
    A_local = params.A_local;     % 3x5 矩阵，动平台铰链点局部坐标
    
    % 初始化力雅可比矩阵
    G_matrix = zeros(6, 6);
    
    %% 1. 计算5个主动力旋量 (Actuation Wrenches)
    for i = 1:5
        % 计算动平台铰链点在全局坐标系下的位置
        A_global = R_mat * A_local(:, i) + P_pos;
        
        % 支链向量与单位方向向量
        leg_vec = A_global - B(:, i);
        leg_len = norm(leg_vec);
        s_i = leg_vec / leg_len; 
        % if i == 3
        %     disp(R_mat)
        %     disp(A_local(:, i))
        % end
        
        % 动平台中心到铰链点的矢量
        r_i = A_global - P_pos;
        
        % 构造单位纯力旋量 (线矢量形式)
        % Plucker坐标: [力方向矢量; 矩矢量 (r x f)]
        moment_i = cross(r_i, s_i);
        G_matrix(:, i) = [s_i; moment_i];
    end
    
    %% 2. 计算SPR支链的约束旋量 (Constraint Wrench)
    % 提取第1条支链(SPR)的单位方向和基座S关节位置
    s_SPR = G_matrix(1:3, 1);       % P关节方向 (支链方向)
    B_SPR = B(:, 1);                % S关节中心在全局系下的坐标 (假设S在基座)
    
    % 获取R关节在全局系下的轴线方向
    u_R_global = params.R_axis_global; 
    
    % 计算约束纯力的方向: f_dir = s x (u_R x s)
    % 该力垂直于支链，且与支链及R轴线共面
    normal_vec = cross(u_R_global, s_SPR);
    f_c_dir = cross(s_SPR, normal_vec);
    
    % 归一化约束力方向
    if norm(f_c_dir) > 1e-6
        f_c_dir = f_c_dir / norm(f_c_dir); 
    else
        % 当 u_R 与 s_SPR 平行时，机构处于极度奇异状态，需异常处理
        warning('SPR支链的R轴与支链方向重合，约束力方向无法唯一确定。');
        f_c_dir = zeros(3,1);
    end
    
    % 构造纯力约束旋量 (Plucker坐标: [力矢量; 动平台中心到力作用线的矩])
    % 因为力穿过基座S关节 (B_SPR)，相对于动平台参考中心 P_pos 的位置矢量为 r_B
    r_B = B_SPR - P_pos; 
    m_c = cross(r_B, f_c_dir); % 约束力对动平台中心产生的力矩
    
    W_constraint = [f_c_dir; m_c];
    
    % 将正确的约束纯力旋量填入雅可比矩阵第6列
    G_matrix(:, 6) = W_constraint;
    
    %% 3. 解算静力平衡方程
    % G * F + W_ext = 0  =>  F = -inv(G) * W_ext
    % 采用左除提高数值稳定性
    if cond(G_matrix) > 1e10
        warning('力雅可比矩阵接近奇异，机构可能处于奇异位形（例如在Home Position发生的活动度问题）。');
    end
    
    Joint_Forces = - G_matrix \ W_ext;

    disp(G_matrix)
    
    % 提取结果
    tau_active = Joint_Forces(1:5);   % 5个主动力
    f_constraint = Joint_Forces(6);   % 约束力/力矩大小
end