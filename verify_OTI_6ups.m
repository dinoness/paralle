%% 计算 6-UPS Stewart 平台在给定参数和位姿下的 OTI
% 参数：r1=50mm, r2=20mm, r3=20*sqrt(2)mm
% 位姿：高度 z=30mm, 姿态 (0°,0°,0°)
% 几何假设：定平台U副均匀分布；动平台S副均匀分布，错开30°，半径60°交替

clear
path_add();
%% 机构参数（单位：mm）
r1 = 50;           % 定平台半径
r2 = 20;           % 动平台内圈半径
r3 = 20*sqrt(2);   % 动平台外圈半径
z  = 30;           % 动平台高度

% 定平台 U 副连接点（6点，均匀分布，0°起始）
alpha = asin(r3/2/r1);
beta = deg2rad(45);
limb_dir = [0-alpha, 0-beta;
            0+alpha, 0+beta;
            deg2rad(120)-alpha, deg2rad(120)-beta;
            deg2rad(120)+alpha, deg2rad(120)+beta;
            deg2rad(240)-alpha, deg2rad(240)-beta;
            deg2rad(240)+alpha, deg2rad(240)+beta;];
B1 = [r1*cos(limb_dir(1,1)); r1*sin(limb_dir(1,1)); 0];
B2 = [r1*cos(limb_dir(2,1)); r1*sin(limb_dir(2,1)); 0];
B3 = [r1*cos(limb_dir(3,1)); r1*sin(limb_dir(3,1)); 0];
B4 = [r1*cos(limb_dir(4,1)); r1*sin(limb_dir(4,1)); 0];
B5 = [r1*cos(limb_dir(5,1)); r1*sin(limb_dir(5,1)); 0];
B6 = [r1*cos(limb_dir(6,1)); r1*sin(limb_dir(6,1)); 0];
B = [B1 B2 B3 B4 B5 B6];

P1_m = [r2*cos(limb_dir(1,2)); r2*sin(limb_dir(1,2)); 0];
P2_m = [r2*cos(limb_dir(2,2)); r2*sin(limb_dir(2,2)); 0];
P3_m = [r2*cos(limb_dir(3,2)); r2*sin(limb_dir(3,2)); 0];
P4_m = [r2*cos(limb_dir(4,2)); r2*sin(limb_dir(4,2)); 0];
P5_m = [r2*cos(limb_dir(5,2)); r2*sin(limb_dir(5,2)); 0];
P6_m = [r2*cos(limb_dir(6,2)); r2*sin(limb_dir(6,2)); 0];
P_m = [P1_m P2_m P3_m P4_m P5_m P6_m];


x_seq = (-50 : 1 : 50);
y_seq = (-50 : 1 : 50);
nx = length(x_seq);
ny = length(y_seq);
m_OTI = zeros(nx, ny);
phi = deg2rad(0);  % 20 
theta = deg2rad(0);  % 10
sigma = deg2rad(0);  % 5

for ix = 1 : nx
    for iy = 1 : ny
        % 姿态 (0°,0°,0°) —— 无旋转，仅 z 方向平移
        R_plant = [cos(phi)*cos(theta)*cos(sigma-phi)-sin(phi)*sin(sigma-phi), -cos(phi)*cos(theta)*sin(sigma-phi)-sin(phi)*cos(sigma-phi), cos(phi)*sin(theta);
                   sin(phi)*cos(theta)*cos(sigma-phi)+cos(phi)*sin(sigma-phi), -sin(phi)*cos(theta)*sin(sigma-phi)+cos(phi)*cos(sigma-phi), sin(phi)*sin(theta);
                   -sin(theta)*cos(sigma-phi),         sin(theta)*sin(sigma-phi),   cos(theta)];
        p_pos = [x_seq(ix); y_seq(iy); z];
        P = zeros(3, 6);
        for i = 1:6
            P(:, i) = R_plant * P_m(:, i) + p_pos;
        end

        %% 计算传递力旋量 ST 与输入运动旋量 SI
        % 6-UPS 支链为 UPS，传递力为沿支链方向的纯力
        ST = zeros(6, 6);
        SI = zeros(6, 6);
        for i = 1:6
            u = P(:,i) - B(:,i);
            u = u / norm(u);  % 支链单位方向向量（由定平台指向动平台）
            % 纯力旋量：力部分 + 力矩部分（参考点为基座标原点）
            ST(:,i) = [u; cross(B(:,i), u)];  % ===========================存疑
            % 输入运动旋量：P 副的纯移动
            SI(:,i) = [zeros(3,1); u];
        end

        %% 计算输出运动旋量 SO
        % 当锁住其余 5 个输入关节时，末端执行器仅剩 1 个自由度。
        % SO(:,i) 应与除第 i 个外的所有传递力旋量互易。
        Omega = [zeros(3,3) eye(3); eye(3) zeros(3,3)];
        SO = zeros(6, 6);
        for i = 1:6
            idx = [1:i-1, i+1:6];
            W = ST(:, idx);  % 6×5 的力旋量矩阵
            SO(:, i) = null(W' * Omega);
            % 调整符号，使第 i 个传递力对输出运动做正功
            if ST(:, i)' * Omega * SO(:, i) < 0
                SO(:, i) = -SO(:, i);
            end
        end

        %% 计算输出传递指标 OTI
        % 依据螺旋理论：eta_i = |$T_i ∘ $O_i| / ||v_p||
        % 其中 v_p 为动平台连接点（力作用点）处的速度
        eta = zeros(6, 1);
        for i = 1:6
            % omega_O = SO(1:3, i);
            % v_O = SO(4:6, i);
            % r_f = P(:, i);  % 力作用点取动平台连接点
            % v_p = v_O + cross(omega_O, r_f);
            % v_p_norm = norm(v_p);
            % if v_p_norm > 1e-10
            %     eta(i) = abs(ST(:, i)' * Omega * SO(:, i)) / v_p_norm;
            % else
            %     eta(i) = 0;
            % end
            eta(i) = screw_efficiency(ST(:, i), SO(:, i), 'p_cstr', P(:, i));
        end

        OTI = min(eta);
        m_OTI(iy, ix) = OTI;
    end
end

%% 输出结果
fprintf('================ 6-UPS OTI 计算结果 ================\n');
fprintf('几何参数：r1=%.2f mm, r2=%.2f mm, r3=%.2f mm, z=%.2f mm\n', r1, r2, r3, z);
fprintf('姿态：(0°, 0°, 0°)\n\n');
fprintf('各支链输出传递系数 eta_i：\n');
for i = 1:6
    fprintf('  eta_%d = %.6f\n', i, eta(i));
end
fprintf('\n输出传递指标 OTI = min(eta_i) = %.6f\n', OTI);
fprintf('====================================================\n');

% figure
fig = figure('Color', [1 1 1]);
% surf(x_seq,y_seq,m_OTI)
% hold on 
contour(x_seq, y_seq, m_OTI);  % 在画图中，x对应的是列，y对应的是行，而在前面的设定中m_OTI(ix,iy)，x是行，y是列，因此出现偏差
xlabel('x')
ylabel('y')
axis equal
