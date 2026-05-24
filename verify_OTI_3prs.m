%% 计算 3-PRS 并联机构在角度空间内的 OTI
% 参数：R1=200mm, R2=130mm, L=220mm
% 角度空间：phi=[-180°, 180°], theta=[0°, 90°]
% 几何假设：定平台与动平台支链 120° 均布，角度分别为 60°, 180°, 300°
% P 副为垂直方向运动，R 副将支链约束在对应竖直平面内

clear
path_add();

%% 机构参数（单位：mm）
R1 = 200;           % 静平台半径
R2 = 130;           % 动平台半径
L  = 220;           % 支链长度

% 支链分布角（定平台与动平台相同）
alpha = deg2rad([60; 180; 300]);

% 定平台 P 副连接点基座标（x,y 固定，z 由运动学确定）
B_xy = [R1*cos(alpha), R1*sin(alpha)];  % 3x2

% 动平台 S 副连接点（动平台坐标系中）
P_m = [R2*cos(alpha), R2*sin(alpha), zeros(3,1)]';  % 3x3

% 竖直平面法向量（R 副轴线方向）
n_plane = [-sin(alpha), cos(alpha), zeros(3,1)]';  % 3x3

%% 角度空间离散化
phi_seq   = deg2rad(-180 : 2 : 180);
theta_seq = deg2rad(0 : 2 : 90);
nphi   = length(phi_seq);
theta_n = length(theta_seq);
m_OTI = zeros(theta_n, nphi);

%% 预分配
eta_all = zeros(3, 1);
eta_saved = false;

for iphi = 1 : nphi
    for itheta = 1 : theta_n
        phi   = phi_seq(iphi);
        theta = theta_seq(itheta);
        
        % T&T 角旋转矩阵 R(phi, theta)
        cp = cos(phi);
        sp = sin(phi);
        ct = cos(theta);
        st = sin(theta);
        
        R = [cp^2*ct + sp^2,     cp*sp*(ct - 1),       cp*st;
             cp*sp*(ct - 1),     sp^2*ct + cp^2,       sp*st;
             -cp*st,             -sp*st,               ct   ];
        
        % 由平面约束（R 副将各支链限制在对应竖直平面内）确定 px, py
        % 推导结果：px = R2*(1-ct)/2 * cos(2*phi)
        %           py = -R2*(1-ct)/2 * sin(2*phi)
        px = R2 * (1 - ct) / 2 * cos(2*phi);
        py = -R2 * (1 - ct) / 2 * sin(2*phi);
        
        % 计算动平台连接点位置（暂令 pz=0，只求水平投影）
        P_temp = R * P_m + [px; py; 0];
        
        % 计算各支链水平投影距离 d_h
        d_h = zeros(3, 1);
        for i = 1:3
            d_h(i) = norm(P_temp(1:2, i) - B_xy(i, :)');
        end
        
        % 若水平投影超过杆长，该位姿不可行
        if any(d_h > L)
            m_OTI(itheta, iphi) = NaN;
            continue;
        end
        
        % 计算动平台 z 坐标 pz，使三个滑块平均高度为 0
        RP_m = R * P_m;
        Pz_m = RP_m(3, :)';  % 3x1
        z_proj = sqrt(L^2 - d_h.^2);  % 上臂构型取正
        pz = mean(z_proj - Pz_m);
        
        % 重新计算动平台连接点
        p_pos = [px; py; pz];
        P = RP_m + p_pos;
        
        % 计算滑块位置（z 坐标）
        B = zeros(3, 3);
        for i = 1:3
            B(:, i) = [B_xy(i, :)'; P(3, i) - z_proj(i)];
        end
        
        %% 计算传递力旋量 ST、输入运动旋量 SI 与平面约束旋量 C
        ST = zeros(6, 3);
        SI = zeros(6, 3);
        C  = zeros(6, 3);
        for i = 1:3
            u = P(:, i) - B(:, i);
            u = u / norm(u);  % 支链单位方向向量（由定平台指向动平台）
            % 纯力旋量：力部分 + 力矩部分（参考点为基座标原点）
            ST(:, i) = [u; cross(B(:, i), u)];
            % 输入运动旋量：P 副为垂直方向纯移动
            SI(:, i) = [zeros(3, 1); 0; 0; 1];
            % 平面约束旋量：R 副将支链限制在竖直平面内
            ni = n_plane(:, i);
            C(:, i) = [ni; cross(P(:, i), ni)];
        end
        
        %% 计算输出运动旋量 SO
        Omega = [zeros(3,3) eye(3); eye(3) zeros(3,3)];
        SO = zeros(6, 3);
        for i = 1:3
            idx = [1:i-1, i+1:3];
            % W 包含：被锁住的 2 条支链的杆长约束 + 所有 3 条支链的平面约束
            W = [ST(:, idx), C];  % 6×5 的约束旋量矩阵
            so_vec = null(W' * Omega);
            if isempty(so_vec) || size(so_vec, 2) < 1
                m_OTI(itheta, iphi) = NaN;
                continue;
            end
            SO(:, i) = so_vec(:, 1);
            % 调整符号，使第 i 个传递力对输出运动做正功
            if ST(:, i)' * Omega * SO(:, i) < 0
                SO(:, i) = -SO(:, i);
            end
        end
        
        %% 计算输出传递指标 OTI
        eta = zeros(3, 1);
        for i = 1:3
            eta(i) = screw_efficiency(ST(:, i), SO(:, i), 'p_cstr', P(:, i));
        end
        
        OTI = min(eta);
        m_OTI(itheta, iphi) = OTI;
        
        % 保存一个中间可行位姿的结果用于输出（phi≈0°, theta≈30°）
        if ~eta_saved && abs(phi) < deg2rad(5) && abs(theta - deg2rad(30)) < deg2rad(5)
            eta_all = eta;
            eta_saved = true;
        end
    end
end

%% 输出结果
fprintf('================ 3-PRS OTI 计算结果 ================\n');
fprintf('几何参数：R1=%.2f mm, R2=%.2f mm, L=%.2f mm\n', R1, R2, L);
fprintf('角度空间：phi=[-180°, 180°], theta=[0°, 90°]\n\n');
fprintf('参考位姿(phi≈0°, theta≈30°)各支链输出传递系数 eta_i：\n');
for i = 1:3
    fprintf('  eta_%d = %.6f\n', i, eta_all(i));
end
fprintf('\n该位姿输出传递指标 OTI = min(eta_i) = %.6f\n', min(eta_all));
fprintf('====================================================\n');

%% 绘图
% fig = figure('Color', [1 1 1]);
% contourf(rad2deg(phi_seq), rad2deg(theta_seq), m_OTI, 20, 'LineColor', 'none');
% hold on
% contour(rad2deg(phi_seq), rad2deg(theta_seq), m_OTI, 5, 'k-', 'LineWidth', 0.5);
% clim([0 1]);
% colorbar;
% colormap(jet);
% xlabel('\phi (°)')
% ylabel('\theta (°)')
% title('3-PRS 输出传递指标 OTI')
% axis tight

fig = figure('Color', [1 1 1]);

%% 1. 构建网格（若 phi_seq/theta_seq 已是网格矩阵，可省略）
[Phi, Theta] = meshgrid(phi_seq, theta_seq);

%% 2. 极坐标 -> 笛卡尔坐标（phi: 方位角[rad], theta: 半径[rad]）
X = Theta .* cos(Phi);
Y = Theta .* sin(Phi);

%% 3. 绘制填充等高线与线框等高线
contourf(X, Y, m_OTI, 20, 'LineColor', 'none');
hold on;
contour(X, Y, m_OTI, 5, 'k-', 'LineWidth', 0.5);

%% 4. （可选）叠加极坐标网格线与刻度，增强可读性
max_r = max(Theta(:));                 % 最大半径（即 theta 最大值）
r_ticks = linspace(0, max_r, 5);       % 径向刻度圈数，可调整
phi_ticks = linspace(-pi, pi, 13);     % 角度线，每 30° 一条

% 绘制同心圆（对应 theta 常数）
for r = r_ticks(2:end)
    th = linspace(0, 2*pi, 200);
    plot(r*cos(th), r*sin(th), '--', 'LineWidth', 0.3, 'Color', [0.5 0.5 0.5]);
end

% 绘制角度射线（对应 phi 常数）
for p = phi_ticks(1:end-1)
    r_vec = linspace(0, max_r, 100);
    plot(r_vec*cos(p), r_vec*sin(p), '--', 'LineWidth', 0.3, 'Color', [0.5 0.5 0.5]);
end

% 径向刻度标签（在最右侧 phi=0° 处标注）
for r = r_ticks(2:end)
    text(r, 0, sprintf(' %.0f°', rad2deg(r)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
        'FontSize', 9, 'Color', [0.2 0.2 0.2]);
end

%% 5. 颜色与格式设置
clim([0 1]);
colorbar;
colormap(jet);
axis equal tight;          % 保证圆形比例正确，消除空白

% 隐藏笛卡尔轴标签（极坐标下不再需要 x/y label）
xlabel('');
ylabel('');
title('3-PRS 输出传递指标 OTI');