clear
path_add();
%% 参数集
%--------parameter3--------
unit_para = 0.001;  % 0.001表示m，1表示mm
basic_paras = basic_read('parameters.xlsx', 'column', 'C');
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


% 局部指数基公式
zeta_r = [0;0;1;0;0;0];  % 旋转基底，在z轴为运动方向的前提下
zeta_p = [0;0;0;0;0;1];  % 平移基底
p_seq = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);
P_v = zeros(3, 5);  % 只变换了方向，没变换起点
P = zeros(3, 5);    % 末端点坐标
% -----end-parameter------


x_seq = (-250 : 10 : 250) * unit_para;
y_seq = (-500 : 10 : 0) * unit_para;
z_seq = -800 * unit_para;
phi_seq = 0;
theta_seq = 0;
Pos_ref_seq = [0*unit_para;0*unit_para;-650*unit_para;0;0];  % line=5 colum=n  角度的单位是° **一列为一组**



% result
nx = length(x_seq);
ny = length(y_seq);
m_ITI = zeros(nx, ny);
m_OTI = zeros(nx, ny);
m_ICI = zeros(nx, ny);
m_OCI = zeros(nx, ny);

for ix = 1 : nx
    for iy = 1 : ny
        Pos_ref = [x_seq(ix); y_seq(iy);z_seq;phi_seq;theta_seq];  % 指定的位姿
        T_ref = pos2trans(Pos_ref(:, 1), B);  % 转移矩阵
        joint_q0 = keni_sol_inverse(T_ref, B, l0_seq, P_m, p_seq);
        P = T_ref(1:3, 4) + T_ref(1:3, 1:3)*P_m;  % T_ref下的各个动平台关节坐标

        % 指标计算
        U = zeros(6, 6, 5);
        ST = zeros(6, 5);  % 传递力旋量
        SC = zeros(6, 1);  % 约束力旋量
        SI = zeros(6, 5);  % 输入运动旋量
        SO = zeros(6, 5);  % 输出运动旋量
        Delta_SO = zeros(6, 1);  % 输出受限运动旋量

        % SPR
        U1 = zeros(6, 6);  % 运动副旋量

        T1 = zeros(4,4,6);
        for i_joint = 1 : 6
            T1(:,:,i_joint) = exp_se3(p_seq(:, i_joint));
        end

        T_temp = eye(4);
        T_zeta = eye(4);
        for i_joint = 1 : 5
            T_temp = T_temp * T_zeta * T1(:,:,i_joint);

            if i_joint ~= 4
                T_zeta = exp_se3(zeta_r * joint_q0(i_joint,1));
            else
                T_zeta = exp_se3(zeta_p * joint_q0(i_joint,1));
            end

            if i_joint ~= 4
                U1(:, i_joint) = adjoint_m(T_temp) * zeta_r;
            else
                U1(:, i_joint) = adjoint_m(T_temp) * zeta_p;
            end
        end
        U(:, :, 1) = U1;
        Omega = [zeros(3,3) eye(3); eye(3) zeros(3,3)];
        SI(:, 1) = U1(:, 4);
        SC(:, 1) = null(U1(:, 1:5)' * Omega);
        WA1 = Omega * U1(:, 1:5) / (U1(:, 1:5)' * U1(:, 1:5));  % 驱动力旋量空间
        TR1 = null(WA1' * Omega);  % 受限运动旋量空间
        ST(:, 1) = [U1(4:6,4); cross(B(:, 1), U1(4:6,4))];  % 纯力，以基座标原点为参考点
        % ST(:, 1) = WA1(:, 4);


        % UPS
        T2 = zeros(4,4,7);
        for i_limb = 2 : 5
            U2 = zeros(6, 6);
            for i_joint = 1 : 7
                T2(:,:,i_joint) = exp_se3(p_seq(:, (7*(i_limb-1) - 1) + i_joint));
            end

            T_temp = eye(4);
            T_zeta = eye(4);
            for i_joint = 1 : 6
                T_temp = T_temp * T_zeta * T2(:,:,i_joint);

                if i_joint ~= 3
                    T_zeta = exp_se3(zeta_r * joint_q0(i_joint,i_limb));
                else
                    T_zeta = exp_se3(zeta_p * joint_q0(i_joint,i_limb));
                end

                if i_joint ~= 3
                    U2(:, i_joint) = adjoint_m(T_temp) * zeta_r;
                else
                    U2(:, i_joint) = adjoint_m(T_temp) * zeta_p;
                end
            end

            U(:, :, i_limb) = U2;
            ST(:, i_limb) = [U2(4:6,3); cross(B(:,i_limb), U2(4:6,3))];
            SI(:, i_limb) = U2(:, 3);
            % ST(:, i_limb) = null([U2(:, 1:2) U2(:, 4:6)]' * Omega);

        end

        %% 计算输出运动旋量 SO，与受限运动旋量Delta_SO
        for i = 1 : 5
            idx = [1:i-1, i+1:5];
            W = [SC, ST(:, idx)];  % 约束力 + 除第i个外的传递力
            SO(:, i) = null(W' * Omega);
            if ST(:, i)' * Omega * SO(:, i) < 0
                SO(:, i) = -SO(:, i);
            end
        end

        Delta_SO(:, 1) = null(ST' * Omega);


        %% ITI OTI
        lambda = zeros(5, 1);
        eta = zeros(5, 1);
        for i = 1: 5
            % ITI
            % lambda(i) = abs(ST(:, i)' * Omega * SI(:, i)) / screw_apparent_power(ST(:, i), SI(:, i), 'p_cstr', P(:, i));  % 分子没有归一化
            lambda(i) = screw_efficiency(ST(:, i), SI(:, i), 'p_cstr', P(:, i));
            % OTI
            % eta(i) = abs(ST(:, i)' * Omega * SO(:, i)) / screw_apparent_power(ST(:, i), SO(:, i), 'p_cstr', P(:, i));  % 分子没有归一化
            eta(i) = screw_efficiency(ST(:, i), SO(:, i), 'p_cstr', P(:, i));
            
        end

        ITI = min(lambda);
        OTI = min(eta);
        LTI = min(ITI, OTI);

        %% ICI OCI
        zeta = screw_efficiency(SC(:, 1), TR1, 'p_cstr', B(:, 1));
        kappa = screw_efficiency(SC(:, 1), Delta_SO(:, 1), 'p_cstr', B(:, 1));
        GCI = min(zeta, kappa);

        m_ITI(iy, ix) = ITI;
        m_OTI(iy, ix) = OTI;
        m_ICI(iy, ix) = zeta;
        m_OCI(iy, ix) = kappa;
    end
end


% figure
fig = figure('Color', [1 1 1]);
s = surf(x_seq,y_seq,m_OTI,'FaceAlpha',0.8);  % 在画图中，x对应的是列，y对应的是行，而在前面的设定中m_OTI(ix,iy)，x是行，y是列，因此出现偏差
xlabel('x')
ylabel('y')
% s.EdgeColor = 'none';
s.FaceColor = 'interp';

[~,~,~, info] = screw_efficiency(SC(:, 1), Delta_SO(:, 1), 'p_cstr', B(:, 1));
