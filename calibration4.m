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

% -----end-struct-parameter------

err = 1e-6;
loop_max = 10;
% ----- input data ------
Pos_m_seq = [0.01;-0.01;-600.02;0.001;-0.001];  % line=5 colum=n
Pos_ref_seq = [0;30;-600;0;10];  % line=5 colum=n  角度的单位是° **一列为一组**
seq_len = length(Pos_ref_seq(1, :));
Pos_err_seq = zeros(5, seq_len);  % 位姿估计误差，优化的目标
Pos_delta_seq = zeros(5, seq_len);  % 位姿扰动序列
% ----- end input data ------

%% 标定步骤
for im = 1 : seq_len
    pos_ref = Pos_ref_seq(:, im);
    % 逆解参考位姿
    joint_q_ref = keni_sol_inverse(pos_ref, limb_dir, B, r1, r2, l0_seq, P_m);

    % 引入真实位姿，判断误差是否小于容差，是则结束，否则进入迭代
    calib_loop = 0;
    while (condition)  % 写判断条件
        if calib_loop > loop_max
            break;
        end
        % 更新运动学参数

        forwrad_keni_loop = 0;
        % 进入正解循环
        while (condition)
            if forwrad_keni_loop > loop_max
                break;
            end

            

            forwrad_keni_loop = forwrad_keni_loop + 1;
        end
        calib_loop = calib_loop + 1;
    end
    
end
