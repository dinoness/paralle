function [T_actual, joint_q_] = keni_sol_forward(joint_q, p_seq, err_max)
% p_seq __ line 螺旋量维数 = 6; colum 各轴排列 = 34 % 1-6,7-13,14-20,21-27,28-34
% joint_q __ line 单轴参数 = 6; colum 总共轴数 = 5 %
if(~exist('err','var'))
    err_max = 1e-6;  % 如果未出现该变量，则对其进行赋值
end

T_actual = zeros(4, 4);
tol = 1e-3;  % 矩阵移项
tol2 = 1e-5;  % 矩阵移项
alpha = 0.5;  % 迭代步长因子
delta_switch = 2.3e-3;
mu = 1e-5; % 初始阻尼因子

joint_q_ = joint_q;
T0 = keni_sol_forward_once(joint_q_, p_seq);
err = err_cal(T0);

J_q = zeros(6,6,5);  % 初始雅可比，其中SPR支链第六列全为0
loop_max = 500;
loop = 0;
J_passive = zeros(6,5,5);  % 去除主动关节后的雅可比矩阵
joint_passive = zeros(24, 1);  % 所有被动关节排列成列向量，也删除了SPR关节的多出的0

err_list = zeros(loop_max,1);
cur_err = norm(err);
while cur_err > err_max  % sum(abs(err))
    % update joint_q
    % 计算雅可比
    J_q = jacobian_space(joint_q_, p_seq);
    
    % KEY_EXPERIENCE 删去了多余的零列
    % 去除主动关节
    J_passive(:,:,1) = [J_q(:,1:3,1) J_q(:,5:6,1) ];
    joint_passive(1: 4) = [joint_q_(1:3, 1); joint_q_(5, 1)];
    for i_limb = 2 : 5
        joint_passive(i_limb*5-5 : i_limb*5-1) = [joint_q_(1:2, i_limb); joint_q_(4:6, i_limb)];
        J_passive(:,:,i_limb) = [J_q(:,1:2, i_limb) J_q(:,4:6,i_limb)];
    end


    J_all = [J_passive(:,1:4,1) -1*J_passive(:,:,2)          zeros(6,5)          zeros(6,5)          zeros(6,5);
             J_passive(:,1:4,1)          zeros(6,5) -1*J_passive(:,:,3)          zeros(6,5)          zeros(6,5);
             J_passive(:,1:4,1)          zeros(6,5)          zeros(6,5) -1*J_passive(:,:,4)          zeros(6,5);
             J_passive(:,1:4,1)          zeros(6,5)          zeros(6,5)          zeros(6,5) -1*J_passive(:,:,5)];

    % 计算步长 (注意这里不再使用固定 alpha)
    delta_q = (J_all'*J_all + mu * eye(24)) \ (J_all' * err);

    max_step = 0.5; % 限制单步最大变化量(弧度或米)，防止飞点
    if norm(delta_q) > max_step
        delta_q = delta_q * (max_step / norm(delta_q));
    end
    joint_passive_new = joint_passive + delta_q;

    % --- 试探性更新 ---
    joint_q_new = joint_q_;
    joint_q_new(1:3, 1) = joint_passive_new(1:3);
    joint_q_new(5, 1) = joint_passive_new(4);
    for i_limb = 2 : 5
        joint_q_new(1:2, i_limb) = joint_passive_new(i_limb*5-5: i_limb*5-4);
        joint_q_new(4:6, i_limb) = joint_passive_new(i_limb*5-3: i_limb*5-1);
    end

    T0_new = keni_sol_forward_once(joint_q_new, p_seq);
    err_new = err_cal(T0_new);

    H = J_all'*J_all;
    g = J_all'*err;
    D = diag(diag(H));
    actual_decrease = norm(cur_err)^2 - norm(err_new)^2;
    expected_decrease = delta_q' * (mu * D * delta_q + g);

    if expected_decrease < 1e-12 
        rho = actual_decrease; 
    else
        rho = actual_decrease / expected_decrease;
    end

    % 判断是否接受该步
    if rho > 0 || norm(err_new) < norm(cur_err)
        joint_passive = joint_passive_new;
        joint_q_ = joint_q_new;
        cur_err = err_new;        
        loop = loop + 1;
        err_list(loop) = norm(cur_err);
        mu = mu * max(1/3, 1 - (2*rho - 1)^3);
    else
        % 误差上升，拒绝该步，增大阻尼以减小步长
        mu = mu * 3;
        % loop 不增加，直接进入下一次循环重新计算 delta_q
        
        % 防止死循环的安全机制
        if mu > 1e5
            warning('阻尼因子过大，可能陷入局部极小值');
            break;
        end
    end



    
    
end

% 输出误差曲线
plot(err_list(1:loop));
if err_list(end) > 2
    disp(err_list(end));
    error("--- 运动学正解未收敛 ---");
end

T_actual = mean(T0, 3);

end


function err = err_cal(T)
    err = zeros(24, 1);
    T_ref = T(:,:,1); % 以支链1的位姿为绝对基准
    for i_limb = 2 : 5
        % 计算其他支链相对于基准支链的偏差
        err(6*(i_limb-2)+1 : 6*(i_limb-2)+6) = log_se3(T(:,:,i_limb) / T_ref);
    end
end