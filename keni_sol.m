% 结构参数形式
p_seq = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);
% 理想结构参数下运动学逆解
joint_q0 = keni_sol_inverse(T_ref, B, l0_seq, P_m, p_seq);
% 结构参数扰动
p_seq2 = p_seq + 0.0001*rand(6,34);
% 求解扰动后的运动学正解
keni_sol_forward(joint_q0, p_seq2, 1e-6)


function T0 = keni_sol_forward(joint_q, p_seq, err_max)
% p_seq __ line 螺旋量维数 = 6; colum 各轴排列 = 34 % 1-6,7-13,14-20,21-27,28-34
% joint_q __ line 单轴参数 = 6; colum 总共轴数 = 5 %
if(~exist('err','var'))
    err_max = 1e-6;  % 如果未出现该变量，则对其进行赋值
end


joint_q_ = joint_q;
T0 = keni_sol_forward_once(joint_q_, p_seq);
err = err_cal(T0);

J_q = zeros(6,6,5);  % 初始雅可比，其中SPR支链第六列全为0
loop_max = 100;
loop = 0;
J_passive = zeros(6,5,5);  % 去除主动关节后的雅可比矩阵
joint_passive = zeros(24, 1);  % 所有被动关节排列成列向量，也删除了SPR关节的多出的0

err_list = zeros(loop_max,1);

while sum(abs(err)) > err_max
    % update joint_q
    % 计算雅可比
    J_q = jacobian_body(joint_q_, p_seq);
    
    % 去除主动关节
    J_passive(:,:,1) = [J_q(:,1:3,1) J_q(:,4:5,1)];
    joint_passive(1: 4) = [joint_q_(1:3, 1); joint_q_(5, 1)];
    for i_limb = 2 : 5
        joint_passive(i_limb*5-5 : i_limb*5-1) = [joint_q_(1:2, i_limb); joint_q_(4:6, i_limb)];
        J_passive(:,:,i_limb) = [J_q(:,1:2, i_limb) J_q(:,4:6,i_limb)];
    end

    J_all = [J_passive(:,1:4,1) -1*J_passive(:,:,2)          zeros(6,5)          zeros(6,5)          zeros(6,5);
                     zeros(6,4)    J_passive(:,:,2) -1*J_passive(:,:,3)          zeros(6,5)          zeros(6,5);
                     zeros(6,4)          zeros(6,5)    J_passive(:,:,3) -1*J_passive(:,:,4)          zeros(6,5);
                     zeros(6,4)          zeros(6,5)          zeros(6,5)    J_passive(:,:,4) -1*J_passive(:,:,5)];
    
    % joint_passive = joint_passive + pinv(J_all'*J_all) * J_all' * err;
    joint_passive = joint_passive + (J_all'*J_all) \ J_all' * err;
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


    
    if(loop < loop_max)
        loop = loop + 1;
    else
        break;
    end

    % fprintf("loop = %d, err = %.4f, rank(J_all'*J_all) = %d\n", loop, sum(abs(err)), rank(J_all'*J_all));
    err_list(loop) = sum(abs(err));
end


plot(err_list);

end



function err = err_cal(T)
    err = zeros(24, 1);
    for i_limb = 1 : 4
        err(6*(i_limb-1)+1 : 6*(i_limb-1)+6) = log_se3(T(:,:,i_limb+1)/T(:,:,i_limb));
    end
end

