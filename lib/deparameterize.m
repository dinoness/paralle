function [B, r1, r2, l0_seq, P_m, limb_dir, joint_u_angle_tilt] = deparameterize(p_seq, P_m_nom, limb_dir_nom)
%DEPARAMETERIZE 从标定后的运动学参数 p_seq 反解几何参数（parameterize 的逆）
%   输入：
%     p_seq        — 6×34 运动学参数矩阵（标定后）
%     P_m_nom      — 3×5 名义动平台铰点坐标，提供 x,y 分量（P_m(1:2,:) 不在 p_seq 中）
%     limb_dir_nom — (可选) 5×1 或 5×2 名义支链方向角，用于填充 SPR 支链（第 1 行）
%   输出：
%     B                  — 3×5 基座铰点坐标
%     r1, r2             — 动平台半径（下/上）
%     l0_seq             — 1×5 初始杆长
%     P_m                — 3×5 动平台铰点坐标（z 分量从 p_seq 更新，x,y 沿用 P_m_nom）
%     limb_dir           — 5×2 支链方向角 [θ_base, θ_move]（rad）；
%                           SPR 支链方向角未编码在 p_seq 中，需由 limb_dir_nom 提供
%     joint_u_angle_tilt — U 副倾角（rad），由 4 条 UPS 支链平均得到

    if nargin < 3
        limb_dir_nom = zeros(5, 2);
    end

    P_m = P_m_nom;
    limb_dir = zeros(5, 2);
    limb_dir(1, :) = limb_dir_nom(1, :);  % SPR 支链方向角沿用名义值

    %% SPR 支链（第 1 列组：p_seq(:, 1:6)）
    T1_1 = exp_se3(p_seq(:, 1));
    B(:, 1) = T1_1(1:3, 4);

    T1_5 = exp_se3(p_seq(:, 5));
    l0_seq(1) = T1_5(3, 4);

    T1_p = exp_se3(p_seq(:, 6));
    P_m(3, 1) = T1_p(1, 4);
    r1_from_spr = T1_p(2, 4);

    %% UPS 支链（第 2–5 列组）
    alpha_vals = zeros(4, 1);
    r1_ups = zeros(2, 1);
    r2_ups = zeros(2, 1);

    for i_limb = 2 : 5
        c0 = 7 * (i_limb - 1);  % 该支链 p_seq 起始列（0-based offset）

        % --- B(:, i_limb)：T01 平移分量 ---
        T01 = exp_se3(p_seq(:, c0 + 0));
        B(:, i_limb) = T01(1:3, 4);

        % --- limb_dir(i_limb, 1) 与 joint_u_angle_tilt：T01 旋转矩阵 z 轴 ---
        % T01(:,3) = [sin(α)·cos(θ); sin(α)·sin(θ); cos(α)]
        z01 = T01(1:3, 3);
        alpha_vals(i_limb - 1) = acos(max(min(z01(3), 1), -1));
        limb_dir(i_limb, 1) = atan2(z01(2), z01(1));

        % --- l0_seq(i_limb)：T34 平移分量 ---
        T34 = exp_se3(p_seq(:, c0 + 3));
        l0_seq(i_limb) = T34(3, 4);

        % --- limb_dir(i_limb, 2)：由 T45 的 dif_limb_dir 反推 ---
        % T45(:,3) = [cos(π+diff); sin(π+diff); 0]
        T45 = exp_se3(p_seq(:, c0 + 4));
        z45 = T45(1:3, 3);
        dif = atan2(z45(2), z45(1)) - pi;
        limb_dir(i_limb, 2) = limb_dir(i_limb, 1) - dif;

        % --- P_m(3,i) 与 r：T_p 平移分量 ---
        % SPR/UPS 支链 2,3: r = r1,  t_p = [P_m(3,i); 0; r1]
        % UPS 支链 4,5:       r = r2,  t_p = [P_m(3,i); 0; r2]
        T_p = exp_se3(p_seq(:, c0 + 6));
        P_m(3, i_limb) = T_p(1, 4);
        if i_limb <= 3
            r1_ups(i_limb - 1) = T_p(3, 4);
        else
            r2_ups(i_limb - 3) = T_p(3, 4);
        end
    end

    joint_u_angle_tilt = mean(alpha_vals);
    r1 = mean([r1_from_spr; r1_ups]);
    r2 = mean(r2_ups);

    % 若名义 limb_dir 为单列（仅 θ_base），则输出也裁剪为单列
    if size(limb_dir_nom, 2) == 1
        limb_dir = limb_dir(:, 1);
    end

end



% function [B, r1, r2, l0_seq, P_m, limb_dir, joint_u_angle_tilt, info] = ...
%     deparameterize(p_seq, P_m_nom, limb_dir_nom, opts_in)
% %DEPARAMETERIZE 从标定后的运动学参数 p_seq 反解几何参数。
% %
% % 该函数是 parameterize 的结构参数投影/反解版本。除保持原接口外，修正了：
% %   1) 用更新后的 r1/r2 与 limb_dir(:,2) 重建 P_m(1:2,:)，避免半径与动平台
% %      铰点坐标不一致；
% %   2) 所有输出方向角均归一化到 [-pi, pi)；
% %   3) 动平台方向角同时利用 T45 与 T_p 中编码的冗余信息；
% %   4) U 副倾角同时利用 T01 与 T12 中编码的冗余信息；
% %   5) 增加尺寸、有限值和半径正定性检查，并返回各支链样本便于诊断。
% %
% % 输入：
% %   p_seq        — 6×34 标定后的运动学参数矩阵
% %   P_m_nom      — 3×5 名义动平台铰点坐标。默认仅用于兼容旧接口；当
% %                  opts.pm_xy_mode = 'nominal' 时用于保留 P_m(1:2,:)
% %   limb_dir_nom — 5×1 或 5×2 名义支链方向角，单位为 rad。SPR 支链的方向角
% %                  未编码在 p_seq 中，因此第 1 行始终沿用该名义值；若只给一列，
% %                  则认为基座侧和动平台侧方向角相同
% %   opts_in      — 可选结构体：
% %                    pm_xy_mode       = 'consistent'（默认）| 'nominal'
% %                                       'consistent'：P_m(1:2,:) 由更新后的
% %                                       r1/r2 和 limb_dir(:,2) 重新生成；
% %                                       'nominal'：保留旧脚本行为，沿用名义 xy
% %                    theta_move_source = 'all'（默认）| 'T45' | 'Tp'
% %                                       'all'：圆平均融合 T45 与 T_p 的估计；
% %                                       'T45'：仅按旧脚本从 T45 提取；
% %                                       'Tp'：仅从 T_p 提取
% %                    alpha_source      = 'all'（默认）| 'T01'
% %                                       'all'：平均 T01 与 T12 的估计；
% %                                       'T01'：仅按旧脚本从 T01 提取
% %
% % 输出：
% %   B                  — 3×5 基座铰点坐标
% %   r1, r2             — 动平台下/上铰点半径
% %   l0_seq             — 1×5 初始杆长
% %   P_m                — 3×5 动平台铰点坐标。默认 xy 与 r1/r2、limb_dir 一致，
% %                        z 从 p_seq 更新
% %   limb_dir           — 5×2 支链方向角 [theta_base, theta_move]，rad，
% %                        已归一化到 [-pi, pi)
% %   joint_u_angle_tilt — 公共 U 副倾角，rad
% %   info               — 投影诊断信息（各支链样本、离散程度及所用选项）
% %
% % 长度单位与 p_seq 保持一致。

%     %% 输入与选项
%     if nargin < 2 || isempty(P_m_nom)
%         P_m_nom = zeros(3, 5);
%     end
%     if nargin < 3 || isempty(limb_dir_nom)
%         limb_dir_nom = zeros(5, 1);
%     end
%     if nargin < 4 || isempty(opts_in)
%         opts_in = struct();
%     end

%     if ~isnumeric(p_seq) || ~isequal(size(p_seq), [6, 34])
%         error('deparameterize:BadPSeqSize', 'p_seq 必须是 6×34 数值矩阵。');
%     end
%     if ~all(isfinite(p_seq(:)))
%         error('deparameterize:NonFinitePSeq', 'p_seq 中存在 NaN 或 Inf。');
%     end
%     if ~isnumeric(P_m_nom) || ~isequal(size(P_m_nom), [3, 5])
%         error('deparameterize:BadPmNomSize', 'P_m_nom 必须是 3×5 数值矩阵。');
%     end
%     if ~isnumeric(limb_dir_nom) || size(limb_dir_nom, 1) ~= 5 || ...
%             ~ismember(size(limb_dir_nom, 2), [1, 2])
%         error('deparameterize:BadLimbDirNomSize', ...
%             'limb_dir_nom 必须是 5×1 或 5×2 数值矩阵。');
%     end
%     if ~isstruct(opts_in)
%         error('deparameterize:BadOptions', 'opts_in 必须是结构体。');
%     end

%     opts = struct( ...
%         'pm_xy_mode', 'consistent', ...
%         'theta_move_source', 'all', ...
%         'alpha_source', 'all');
%     opt_names = fieldnames(opts_in);
%     for k = 1:numel(opt_names)
%         if ~isfield(opts, opt_names{k})
%             error('deparameterize:UnknownOption', '未知选项：%s。', opt_names{k});
%         end
%         opts.(opt_names{k}) = opts_in.(opt_names{k});
%     end
%     opts.pm_xy_mode = lower(char(opts.pm_xy_mode));
%     opts.theta_move_source = lower(char(opts.theta_move_source));
%     opts.alpha_source = lower(char(opts.alpha_source));

%     if ~ismember(opts.pm_xy_mode, {'consistent', 'nominal'})
%         error('deparameterize:BadPmXYMode', ...
%             'opts.pm_xy_mode 只能为 ''consistent'' 或 ''nominal''。');
%     end
%     if ~ismember(opts.theta_move_source, {'all', 't45', 'tp'})
%         error('deparameterize:BadThetaMoveSource', ...
%             'opts.theta_move_source 只能为 ''all''、''T45'' 或 ''Tp''。');
%     end
%     if ~ismember(opts.alpha_source, {'all', 't01'})
%         error('deparameterize:BadAlphaSource', ...
%             'opts.alpha_source 只能为 ''all'' 或 ''T01''。');
%     end

%     %% 初始化
%     B = zeros(3, 5);
%     P_m = zeros(3, 5);
%     l0_seq = zeros(1, 5);
%     limb_dir = zeros(5, 2);

%     limb_dir_nom = wrap_to_pi(limb_dir_nom);
%     if size(limb_dir_nom, 2) == 1
%         limb_dir_nom = [limb_dir_nom, limb_dir_nom];
%     end
%     limb_dir(1, :) = limb_dir_nom(1, :);  % SPR 支链方向角未编码在 p_seq 中

%     alpha_limb = zeros(4, 1);
%     alpha_sources = zeros(4, 2);       % [T01, T12]
%     theta_move_sources = zeros(4, 3);  % [T45, Tp-x列, Tp-y列]
%     r1_samples = zeros(3, 1);          % SPR + UPS 2,3
%     r2_samples = zeros(2, 1);          % UPS 4,5

%     %% SPR 支链：p_seq(:, 1:6)
%     T1_1 = exp_se3(p_seq(:, 1));
%     B(:, 1) = T1_1(1:3, 4);

%     T1_5 = exp_se3(p_seq(:, 5));
%     l0_seq(1) = T1_5(3, 4);

%     T1_p = exp_se3(p_seq(:, 6));       % t = [P_m(3,1); r1; 0]
%     P_m(3, 1) = T1_p(1, 4);
%     r1_samples(1) = T1_p(2, 4);

%     %% UPS 支链：p_seq(:, 7:34)
%     for i_limb = 2:5
%         c0 = 7 * (i_limb - 1);          % 当前支链在 p_seq 中的起始列
%         k = i_limb - 1;                 % 1..4，用于 UPS 样本索引

%         T01 = exp_se3(p_seq(:, c0 + 0));
%         T12 = exp_se3(p_seq(:, c0 + 1));
%         T34 = exp_se3(p_seq(:, c0 + 3));
%         T45 = exp_se3(p_seq(:, c0 + 4));
%         T_p = exp_se3(p_seq(:, c0 + 6));

%         % 基座铰点与初始杆长
%         B(:, i_limb) = T01(1:3, 4);
%         l0_seq(i_limb) = T34(3, 4);

%         % T01(:,3) = [sin(alpha)*cos(theta_b);
%         %             sin(alpha)*sin(theta_b);
%         %             cos(alpha)]
%         z01 = T01(1:3, 3);
%         alpha_from_T01 = atan2(hypot(z01(1), z01(2)), z01(3));
%         theta_base = atan2(z01(2), z01(1));

%         % T12(:,1) = [0; sin(pi-alpha); cos(pi-alpha)]，因此
%         % alpha = atan2(T12(2,1), -T12(3,1))
%         x12 = T12(1:3, 1);
%         alpha_from_T12 = atan2(x12(2), -x12(3));
%         alpha_from_T12 = min(max(alpha_from_T12, 0), pi);

%         alpha_sources(k, :) = [alpha_from_T01, alpha_from_T12];
%         switch opts.alpha_source
%             case 'all'
%                 alpha_limb(k) = 0.5 * (alpha_from_T01 + alpha_from_T12);
%             case 't01'
%                 alpha_limb(k) = alpha_from_T01;
%         end

%         limb_dir(i_limb, 1) = wrap_to_pi(theta_base);

%         % T45(:,3) = [cos(pi+dif); sin(pi+dif); 0]
%         % dif = theta_base - theta_move
%         z45 = T45(1:3, 3);
%         dif_limb_dir = wrap_to_pi(atan2(z45(2), z45(1)) - pi);
%         theta_move_T45 = wrap_to_pi(theta_base - dif_limb_dir);

%         % T_p 的旋转矩阵也直接编码 theta_move：
%         % T_p(:,1) = [0; sin(theta_m); -cos(theta_m)]
%         % T_p(:,2) = [0; -cos(theta_m); -sin(theta_m)]
%         R_p = T_p(1:3, 1:3);
%         theta_move_Tp_x = atan2(R_p(2, 1), -R_p(3, 1));
%         theta_move_Tp_y = atan2(-R_p(3, 2), -R_p(2, 2));

%         theta_move_sources(k, :) = ...
%             wrap_to_pi([theta_move_T45, theta_move_Tp_x, theta_move_Tp_y]);
%         switch opts.theta_move_source
%             case 'all'
%                 theta_move = circ_mean(theta_move_sources(k, :));
%             case 't45'
%                 theta_move = theta_move_T45;
%             case 'tp'
%                 theta_move = circ_mean([theta_move_Tp_x, theta_move_Tp_y]);
%         end
%         limb_dir(i_limb, 2) = wrap_to_pi(theta_move);

%         % T_p 平移：t = [P_m(3,i); 0; r]
%         P_m(3, i_limb) = T_p(1, 4);
%         if i_limb <= 3
%             r1_samples(i_limb) = T_p(3, 4);       % i_limb=2,3 -> 样本 2,3
%         else
%             r2_samples(i_limb - 3) = T_p(3, 4);   % i_limb=4,5 -> 样本 1,2
%         end
%     end

%     %% 将各支链样本投影回 parameterize 所需的公共结构参数
%     joint_u_angle_tilt = mean(alpha_limb);
%     r1 = mean(r1_samples);
%     r2 = mean(r2_samples);

%     if r1 <= 0 || r2 <= 0
%         error('deparameterize:NonPositiveRadius', ...
%             ['反解得到的平台半径非正：r1=%.12g, r2=%.12g。', ...
%              '请检查 p_seq 或标定结果。'], r1, r2);
%     end
%     if any(l0_seq <= 0)
%         warning('deparameterize:NonPositiveL0', ...
%             '存在非正初始杆长，请检查 p_seq 或标定结果。');
%     end

%     %% 重建动平台铰点坐标
%     % 推荐默认：P_m(1:2,:) 与更新后的 r1/r2、limb_dir(:,2) 严格一致。
%     % 若仅为复现旧脚本输出，可设置 opts.pm_xy_mode='nominal'。
%     for i_limb = 1:5
%         if i_limb <= 3
%             r_i = r1;
%         else
%             r_i = r2;
%         end

%         switch opts.pm_xy_mode
%             case 'consistent'
%                 P_m(1, i_limb) = r_i * cos(limb_dir(i_limb, 2));
%                 P_m(2, i_limb) = r_i * sin(limb_dir(i_limb, 2));
%             case 'nominal'
%                 P_m(1:2, i_limb) = P_m_nom(1:2, i_limb);
%         end
%     end

%     limb_dir = wrap_to_pi(limb_dir);

%     %% 诊断信息
%     info = struct();
%     info.options = opts;
%     info.alpha_limb = alpha_limb;
%     info.alpha_sources = alpha_sources;
%     info.r1_samples = r1_samples;
%     info.r2_samples = r2_samples;
%     info.theta_move_sources = theta_move_sources;
%     info.alpha_spread = max(alpha_limb) - min(alpha_limb);
%     info.r1_spread = max(r1_samples) - min(r1_samples);
%     info.r2_spread = max(r2_samples) - min(r2_samples);
%     dtheta = wrap_to_pi(theta_move_sources - theta_move_sources(:, 1));
%     info.theta_move_spread = max(abs(dtheta), [], 2);
% end


% function angle_out = wrap_to_pi(angle_in)
% %WRAP_TO_PI 将角度归一化到 [-pi, pi)。
%     angle_out = mod(angle_in + pi, 2*pi) - pi;
% end


% function angle_mean = circ_mean(angle_seq)
% %CIRC_MEAN 圆平均，避免直接平均导致的跨 ±pi 错误。
%     angle_seq = angle_seq(:).';
%     s = mean(sin(angle_seq));
%     c = mean(cos(angle_seq));
%     if hypot(s, c) < 1e-12
%         % 角度几乎完全反向时圆平均不唯一，退回第一项并保留诊断样本。
%         angle_mean = angle_seq(1);
%     else
%         angle_mean = atan2(s, c);
%     end
%     angle_mean = wrap_to_pi(angle_mean);
% end