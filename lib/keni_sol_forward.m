function T0 = keni_sol_forward(joint_q, p_seq, err_max)
% p_seq __ line 螺旋量维数 = 6; colum 各轴排列 = 34 % 1-6,7-13,14-20,21-27,28-34
% joint_q __ line 单轴参数 = 6; colum 总共轴数 = 5 %
if(~exist('err','var'))
    err_max = 1e-6;  % 如果未出现该变量，则对其进行赋值
end


joint_q_ = joint_q;
T0 = keni_sol_forward_once(joint_q_, p_seq);
err = err_cal(T0);
J_q = zeros(6,6,5);

loop_max = 100;
loop = 0;

while sum(abs(err)) > err_max
    % update joint_q
    % 计算雅可比
    J_q = jacobian(joint_q_, p_seq);
    joint_passive = zeros(25, 1);
    J_passive = zeros(6,5,5);


    joint_passive(1: 5) = [joint_q_(1:3, 1); joint_q_(5:6, 1)];
    J_passive(:,:,1) = [J_q(:,1:3,1) J_q(:,5:6,1)];
    for i_limb = 2 : 5
        joint_passive(i_limb*5-4 : i_limb*5) = [joint_q_(1:2, i_limb); joint_q_(4:6, i_limb)];
        J_passive(:,:,i_limb) = [J_q(:,1:2, i_limb) J_q(:,4:6,i_limb)];
    end

    % 第5列为全零，导致奇异
    J_all = [J_passive(:,:,1) -1*J_passive(:,:,2)          zeros(6,5)          zeros(6,5)          zeros(6,5);
                   zeros(6,5)    J_passive(:,:,2) -1*J_passive(:,:,3)          zeros(6,5)          zeros(6,5);
                   zeros(6,5)          zeros(6,5)    J_passive(:,:,3) -1*J_passive(:,:,4)          zeros(6,5);
                   zeros(6,5)          zeros(6,5)          zeros(6,5)    J_passive(:,:,4) -1*J_passive(:,:,5)];
    
    joint_passive = joint_passive + (J_all'*J_all) \ J_all' * err;
    disp("joint_passive");
    disp(joint_passive);

    rank(J_all'*J_all)


    joint_q_(1:3, 1) = joint_passive(1:3);
    joint_q_(5:6, 1) = joint_passive(4:5);
    for i_limb = 2 : 5
        joint_q_(1:2, i_limb) = joint_passive(i_limb*5-4: i_limb*5-3);
        joint_q_(4:6, i_limb) = joint_passive(i_limb*5-2: i_limb*5);
    end


    T0 = keni_sol_forward_once(joint_q_, p_seq);
    err = err_cal(T0);


    % disp("err");
    % disp(err);
    if(loop < loop_max)
        loop = loop + 1;
    else
        break;
    end
end


end



function err = err_cal(T)
    err = zeros(24, 1);
    for i_limb = 1 : 4
        err(6*(i_limb-1)+1 : 6*(i_limb-1)+6) = log_se3(T(:,:,i_limb+1)/T(:,:,i_limb));
    end
end

