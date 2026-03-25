function [T_actual, joint_q_] = keni_sol_forward(joint_q, p_seq, err_max)
% p_seq __ line 螺旋量维数 = 6; colum 各轴排列 = 34 % 1-6,7-13,14-20,21-27,28-34
% joint_q __ line 单轴参数 = 6; colum 总共轴数 = 5 %
if(~exist('err','var'))
    err_max = 1e-6;  % 如果未出现该变量，则对其进行赋值
end

T_actual = zeros(4, 4);
tol = 1e-5;  % 矩阵移项
alpha = 0.5;  % 迭代步长因子


joint_q_ = joint_q;
T0 = keni_sol_forward_once(joint_q_, p_seq);
err = err_cal(T0);

J_q = zeros(6,6,5);  % 初始雅可比，其中SPR支链第六列全为0
loop_max = 50;
loop = 0;
J_passive = zeros(6,5,5);  % 去除主动关节后的雅可比矩阵
joint_passive = zeros(24, 1);  % 所有被动关节排列成列向量，也删除了SPR关节的多出的0

err_list = zeros(loop_max,1);

while norm(err) > err_max  % sum(abs(err))
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

    % if loop == 0
    %     % disp(J_q)
    %     % disp(J_passive)
    %     % disp(joint_q_)
    %     % disp(joint_passive)
    %     disp(err)
    %     for i_limb = 1 : 4
    %         exp_se3(err(6*(i_limb-1)+1 : 6*(i_limb-1)+6))
    %     end

    % end

    J_all = [J_passive(:,1:4,1) -1*J_passive(:,:,2)          zeros(6,5)          zeros(6,5)          zeros(6,5);
                     zeros(6,4)    J_passive(:,:,2) -1*J_passive(:,:,3)          zeros(6,5)          zeros(6,5);
                     zeros(6,4)          zeros(6,5)    J_passive(:,:,3) -1*J_passive(:,:,4)          zeros(6,5);
                     zeros(6,4)          zeros(6,5)          zeros(6,5)    J_passive(:,:,4) -1*J_passive(:,:,5)];
    
    % joint_passive = joint_passive + alpha * pinv(J_all'*J_all) * J_all' * err;
    joint_passive = joint_passive + alpha * ((J_all'*J_all + tol*eye(24)) \ (J_all' * err));
    % joint_passive = joint_passive + alpha * (J_all'*J_all) \ J_all' * err;
    % joint_passive = joint_passive + J_all \ err;
    % disp("joint_passive");
    % disp(joint_passive);

    % disp("rank J_all'*J_all");
    % disp(rank(J_all'*J_all));


    % 还原为矩阵模式
    joint_q_(1:3, 1) = joint_passive(1:3);
    joint_q_(5, 1) = joint_passive(4);
    for i_limb = 2 : 5
        joint_q_(1:2, i_limb) = joint_passive(i_limb*5-5: i_limb*5-4);
        joint_q_(4:6, i_limb) = joint_passive(i_limb*5-3: i_limb*5-1);
    end


    T0 = keni_sol_forward_once(joint_q_, p_seq);
    err = err_cal(T0);
    

    % 记录循环次数
    if(loop < loop_max)
        loop = loop + 1;
        err_list(loop) = norm(err);  % sum(abs(err))
    else
        break;
    end

    % fprintf("loop = %d, err = %.4f, rank(J_all'*J_all) = %d\n", loop, sum(abs(err)), rank(J_all'*J_all));
    
    
end

% 输出误差曲线
% plot(err_list(1:loop));
if err_list(end) > 2
    disp(err_list(end));
    error("--- 运动学正解未收敛 ---");
end

T_actual = mean(T0, 3);

end



function err = err_cal(T)
    err = zeros(24, 1);
    for i_limb = 1 : 4
        if rcond(T(:,:,i_limb)) < 1e-14
            % disp(T(:,:,i_limb));
        end
        err(6*(i_limb-1)+1 : 6*(i_limb-1)+6) = log_se3(T(:,:,i_limb+1)/T(:,:,i_limb));
    end
end

