% compute_oti.m — 输入结构参数和末端位姿，计算该位姿下的 OTI（输出传递指标）
%
% 使用方式：
%   1. 直接运行（使用内置示例参数）
%   2. 修改下方 USER INPUT 区域的参数后运行
%
% 输出指标说明：
%   OTI — 输出传递指标 (Output Transmission Index)，OTI = min(eta_i)
%         eta_i 为第 i 条支链的传递力旋量与输出运动旋量之间的螺旋效率
%         OTI ∈ [0,1]，越接近 1 表示力/运动传递性能越好
%   ITI — 输入传递指标 (Input Transmission Index)，ITI = min(lambda_i)
%   LTI — 局部传递指标 (Local Transmission Index)，LTI = min(ITI, OTI)
%   GCI — 全局约束指标 (Global Constraint Index)
%
% 参考文献：c02 并联机构的运动学误差建模及参数可辨识性分析_孔令雨.pdf

clear
path_add();
fprintf('>>>= OTI 计算开始 (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));

%% ==================== USER INPUT ====================
% 选择参数来源：'excel' 从 parameters.xlsx 读取，'manual' 手动输入
param_source = 'excel';

% ── 末端位姿 ──
% 格式：[x; y; z; phi; theta]
%   x, y, z    — 动平台中心位置 (m)
%   phi        — z轴在xOy平面投影与x轴的夹角 (rad, 当 unit='rad')
%   theta      — z轴与世界z轴 [0;0;1] 的夹角 (rad, 当 unit='rad')
% 注意：下方 unit 参数决定 phi/theta 的单位是 'rad' 还是 'deg'
pose = [0; 0.1; -0.9; 0; 0];          % 位姿向量 (SI: m, rad)
pose_unit = 'rad';                     % 位姿角度单位: 'rad' 或 'deg'

%% ==================== 参数定义 ====================

if strcmp(param_source, 'excel')
    % 从 parameters.xlsx 读取结构参数
    basic_paras = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'm');
    unit_para = basic_paras.unit_para;

    r1 = basic_paras.r1;
    r2 = basic_paras.r2;
    limb_dir = basic_paras.limb_dir;
    B = basic_paras.B;
    P_m = basic_paras.P_m;
    l0_seq = basic_paras.l0_seq;
    joint_u_angle_tilt = basic_paras.joint_u_angle_tilt;

    fprintf('\n── 结构参数来源: parameters.xlsx (column B) ──\n');

else
    % 手动输入结构参数（全部使用 SI 单位: m, rad）
    unit_para = 0.001;  % 1 表示 mm 输入，0.001 表示 m 输入

    % 基座参数
    R1 = 550 * unit_para;    % 基座下半径
    R2 = 500 * unit_para;    % 基座上半径
    H  = 0   * unit_para;    % 基座上下层高差

    % 动平台参数
    r1 = 100 * unit_para;    % 动平台下半径
    r2 = 80  * unit_para;    % 动平台上半径
    h  = 10  * unit_para;    % 动平台上下层高差
    L_tool = 0 * unit_para;  % 工具长度

    % 支链方向角 [θ_base, θ_move] (rad)
    % 列1: 基座侧方位角, 列2: 动平台侧方位角
    joint_u_angle_tilt = deg2rad(5);  % U 副倾角
    limb_dir = deg2rad([...
        90,  90;     % 支链1 (SPR)
        210, 210;    % 支链2 (UPS)
        330, 330;    % 支链3 (UPS)
        30,  30;     % 支链4 (UPS)
        150, 150;    % 支链5 (UPS)
    ]);

    % 初始杆长 (m)
    l0_seq = [800; 800; 800; 800; 800] * unit_para;

    % 计算基座铰点 B
    B1 = [R1*cos(limb_dir(1,1)); R1*sin(limb_dir(1,1)); 0];
    B2 = [R1*cos(limb_dir(2,1)); R1*sin(limb_dir(2,1)); 0];
    B3 = [R1*cos(limb_dir(3,1)); R1*sin(limb_dir(3,1)); 0];
    B4 = [R2*cos(limb_dir(4,1)); R2*sin(limb_dir(4,1)); H];
    B5 = [R2*cos(limb_dir(5,1)); R2*sin(limb_dir(5,1)); H];
    B = [B1 B2 B3 B4 B5];

    % 计算动平台铰点 P_m
    P1_m = [r1*cos(limb_dir(1,2)); r1*sin(limb_dir(1,2)); L_tool];
    P2_m = [r1*cos(limb_dir(2,2)); r1*sin(limb_dir(2,2)); L_tool];
    P3_m = [r1*cos(limb_dir(3,2)); r1*sin(limb_dir(3,2)); L_tool];
    P4_m = [r2*cos(limb_dir(4,2)); r2*sin(limb_dir(4,2)); L_tool+h];
    P5_m = [r2*cos(limb_dir(5,2)); r2*sin(limb_dir(5,2)); L_tool+h];
    P_m = [P1_m P2_m P3_m P4_m P5_m];

    fprintf('\n── 结构参数来源: 手动输入 ──\n');
end

%% ==================== 位姿转换 ====================
% pos2trans 将 [x;y;z;phi;theta] 转换为 4×4 齐次变换矩阵
T_ref = pos2trans(pose, B, 'unit', pose_unit);

fprintf('\n── 末端位姿 ──\n');
fprintf('  位置: [%.4f, %.4f, %.4f] mm\n', pose(1)/unit_para, pose(2)/unit_para, pose(3)/unit_para);
if strcmp(pose_unit, 'deg')
    fprintf('  姿态: phi=%.4f°, theta=%.4f°\n', pose(4), pose(5));
else
    fprintf('  姿态: phi=%.4f°, theta=%.4f°\n', rad2deg(pose(4)), rad2deg(pose(5)));
end

%% ==================== 运动学参数化 ====================
p_seq = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

%% ==================== 计算 OTI ====================
[LTI, GCI, OTI, ITI, info] = compute_ltigci(T_ref, B, l0_seq, P_m, p_seq);

%% ==================== 结果输出 ====================
fprintf('\n========== OTI 计算结果 ==========\n');
fprintf('  OTI (输出传递指标) = %.6f  (min eta_i)\n', OTI);
fprintf('  ITI (输入传递指标) = %.6f  (min lambda_i)\n', ITI);
fprintf('  LTI (局部传递指标) = %.6f  (min(ITI, OTI))\n', LTI);
fprintf('  GCI (全局约束指标) = %.6f\n', GCI);
fprintf('===================================\n');

fprintf('\n>>>= OTI 计算完成 (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));
