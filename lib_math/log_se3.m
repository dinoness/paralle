function xi = log_se3(T)
    % SE3_to_se3: 将 4x4 变换矩阵 T 映射为 6x1 李代数向量 xi = [rho; omega]
    % 输入: T - 4x4 齐次变换矩阵
    % 输出: xi - 6x1 李代数向量 [rho1, rho2, rho3, omega1, omega2, omega3]'

    R = T(1:3, 1:3);
    t = T(1:3, 4);

    % 1. 计算旋转角 theta
    tr_R = trace(R);
    theta = acos((tr_R - 1) / 2);

    % 2. 处理特殊情况：旋转极小时 (theta -> 0)
    if abs(theta) < 1e-6
        omega = [0; 0; 0];
        V_inv = eye(3);
    elseif abs(theta - pi) < 1e-6
        % 情况 2: 旋转接近 180度 (theta -> pi)
        % 此时 R = 2*v*v' - I, 提取旋转轴 v
        v2 = [ (R(1,1)+1)/2; (R(2,2)+1)/2; (R(3,3)+1)/2 ];
        v2(v2 < 0) = 0; % 防止数值误差导致负数
        v = sqrt(v2);
        
        % 确定符号
        if v(1) > 0
            v(2) = v(2) * sign(R(1,2));
            v(3) = v(3) * sign(R(1,3));
        elseif v(2) > 0
            v(3) = v(3) * sign(R(2,3));
        end
        
        omega = theta * v;
        
        % 计算 V_inv
        omega_skew = [0, -omega(3), omega(2);
                      omega(3), 0, -omega(1);
                     -omega(2), omega(1), 0];
        V_inv = eye(3) - 0.5 * omega_skew + ...
                (1/theta^2 * (1 - (theta * sin(theta)) / (2 * (1 - cos(theta))))) * (omega_skew^2);
    else
        % 计算旋转向量的方向单位向量 (利用 R - R')
        omega_hat_skew = (theta / (2 * sin(theta))) * (R - R');
        omega = [omega_hat_skew(3, 2); omega_hat_skew(1, 3); omega_hat_skew(2, 1)];
        
        % 3. 计算 V 的逆矩阵，用于求解 rho
        % V = I + (1-cos(theta))/theta^2 * omega_hat + (theta-sin(theta))/theta^3 * omega_hat^2
        % 这里直接给出 V_inv 的闭式解以提高效率
        omega_skew = [0, -omega(3), omega(2);
                      omega(3), 0, -omega(1);
                     -omega(2), omega(1), 0];
        
        V_inv = eye(3) - 0.5 * omega_skew + ...
                (1/theta^2 * (1 - (theta * sin(theta)) / (2 * (1 - cos(theta))))) * (omega_skew^2);
    end

    % 4. 计算 rho
    rho = V_inv * t;

    % 5. 组合成 6 维向量
    xi = [omega; rho];
end