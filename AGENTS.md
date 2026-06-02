# AGENTS.md — SPR-4UPS Parallel Mechanism Kinematic Analysis

## Project Overview

This is a MATLAB project for kinematic modeling, calibration, workspace analysis, force/velocity analysis, and structural parameter optimization of a **SPR-4UPS parallel mechanism** (5-DOF). The mechanism consists of one SPR limb (Spherical-Prismatic-Revolute) and four UPS limbs (Universal-Prismatic-Spherical).

All kinematics are formulated using **screw theory** and the **Product of Exponentials (POE)** formula on Lie groups SE(3) / se(3). The primary reference is the included PDF (孔令雨, "并联机构的运动学误差建模及参数可辨识性分析").

**Code language:** MATLAB (all comments and documentation in Chinese).

## Repository Structure

```
paralle_/
├── calibration2/3/4.m        # 运动学标定 (LM 阻尼最小二乘, 第3章)
├── calibration5.m            # 运动学标定 (总体最小二乘法 TLS, 第4章)
├── parameter_optimize.m      # 单层参数扫描与性能指标作图
├── parameter_optimize2.m     # 结构参数优化 (parameter optimization)
├── parameter_optimize_surrogate.m  # 基于 surrogateopt 的参数优化
├── workspace_discrete.m/v2/v3.m    # 圆柱工作空间离散搜索 (cylindrical workspace search with OTI/GCI)
├── Joint_angle_search.m      # 关节角搜索与工作空间可视化
├── force_analysis_screw.m    # 旋量理论力分析 (screw theory based force analysis)
├── velocity_static_force.m   # 速度与静力分析
├── structure_err.m           # 结构误差对位姿影响分析与可视化
├── scerw_theory.m            # 旋量理论验证
├── screw_verify.m            # 旋量计算结果验证
├── keni_sol.m                # 运动学正/逆解验证
├── controller_verify.m       # 控制器验证
├── verify_OTI_3prs.m         # 3PRS 机构 OTI 验证
├── verify_OTI_6ups.m         # 6UPS 机构 OTI 验证
├── test.m / test_script.m    # 测试脚本
├── path_add.m                # 路径初始化 (adds all lib dirs to MATLAB path)
├── parameters.xlsx           # 机构几何参数表 (多列对应不同参数集)
├── optimization_result.mat   # 优化结果缓存
├── note.md                   # MATLAB 语法备忘
│
├── lib/                      # 核心运动学与标定库
│   ├── parameterize.m        # 几何参数 → 旋量序列 (p_seq, xi_seq)
│   ├── keni_sol_forward.m    # 运动学正解 (Newton-Raphson 迭代)
│   ├── keni_sol_forward_once.m  # 单步正运动学 (POE 正解)
│   ├── keni_sol_inverse.m    # 运动学逆解 (Paden-Kahan 子问题)
│   ├── jacobian_body.m       # 物体雅可比
│   ├── jacobian_space.m      # 空间雅可比
│   ├── pos2trans.m / trans2pos.m  # 位姿向量 ↔ 齐次变换矩阵
│   ├── spr4ups_build_limb_twists.m  # 构建支链关节运动旋量
│   ├── build_line_wrench.m   # 构建线力旋量
│   ├── svd_nullspace.m       # SVD 零空间
│   ├── null_rowspace_z.m     # 行空间/零空间分解
│   ├── transform_matrix_cal.m  # 变换矩阵计算
│   ├── calib_iter_matrix.m / calib_iter_matrix2.m  # 标定迭代矩阵 Jp_bar
│   ├── calib_iter_row_space_matrix.m  # 行空间分解 (U, N, V_prep, M)
│   ├── calib_iter_restore_matrix.m    # 恢复矩阵 (Lambda, Ap)
│   ├── restore_full_param_increment.m # 参数增量恢复为完整形式
│   └── rebuild_xi_seq_from_p.m        # p_seq 重建 xi_seq
│
├── lib_math/                 # 数学基础库
│   ├── exp_se3.m / log_se3.m    # SE(3) 指数/对数映射
│   ├── exp_so3.m / log_so3.m    # SO(3) 指数/对数映射
│   ├── adjoint_m.m              # 伴随变换
│   ├── skew.m                   # 反对称矩阵
│   ├── Paden_Kahan1.m           # PK 子问题 1 (绕单轴旋转)
│   ├── Paden_Kahan2.m           # PK 子问题 2 (绕两相交轴旋转)
│   ├── recip.m                  # 互易积
│   ├── screw_unitize.m          # 旋量单位化
│   ├── screw_apparent_power.m   # 旋量视在功率
│   └── screw_efficiency.m       # 旋量效率
│
├── lib_para/                 # 参数读写
│   └── basic_read.m          # 从 Excel 读取机构参数，支持列选择与单位设置
│
├── lib_opt_struct/           # 结构优化
│   ├── compute_ltigci.m      # 计算 LTI, GCI, OTI, ITI 性能指标
│   └── evaluate_spr4ups_objective.m  # surrogateopt 目标函数
│
├── lib_calib/                # 标定辅助库
│   └── solve_tls.m           # 总体最小二乘法 (TLS) 求解器，增广矩阵 SVD 方法
│
├── AI_code/                  # AI 辅助生成代码
│   └── SPR_4UPS_Statics.m    # SPR-4UPS 静力分析函数
│
├── other/                    # 工具/探索/旧版代码
│   ├── Workspace.m / Workspace2.m   # 旧版工作空间分析
│   ├── Joint_angle_search_1016.m    # 旧版关节搜索
│   ├── kinematics_test.m / base_verify.m / view_test.m  # 验证测试
│   ├── Circle.m / ArcDraw.m         # 绘图辅助
│   ├── frame_assistant_cal.m        # 坐标系辅助计算
│   └── t1.m / t2.m                  # 临时测试
│
└── sim_para_process/         # RecurDyn 仿真数据处理
    ├── recurdyn0427.m
    ├── with_gravity.csv
    └── no_gravity.csv
```

## No Build System

There is **no build system or package manager**. This is pure MATLAB — each `.m` file is a script or function. To run any analysis, open MATLAB, set the working directory to the project root, and run the desired script. Most scripts call `path_add()` at the top to set up the MATLAB search path, or call `addpath(genpath('./lib'))` directly.

## How to Run

1. Open MATLAB
2. Set current folder to the project root (`paralle_/`)
3. Run one of the top-level scripts, e.g.:
   - `>> workspace_discrete_v3` — workspace analysis
   - `>> calibration4` — kinematic calibration (Levenberg-Marquardt)
   - `>> calibration5` — kinematic calibration (Total Least Squares, Chapter 4)
   - `>> parameter_optimize_surrogate` — structural optimization
   - `>> structure_err` — structural error visualization
4. Most scripts call `path_add()` (or equivalent `addpath`) at the top

## Parameter Management

All mechanism geometry is stored in `parameters.xlsx`. The file has **multiple columns** (B, C, …, M), each representing a different parameter set. Columns encode:

| Row | Parameter | Unit in Excel |
|-----|-----------|--------------|
| 1 | l_max | mm |
| 2 | l_min | mm |
| 3 | R1 (base radius, lower platform) | mm |
| 4 | R2 (base radius, upper platform) | mm |
| 5 | H (base vertical offset) | mm |
| 6 | r1 (moving platform radius, lower) | mm |
| 7 | r2 (moving platform radius, upper) | mm |
| 8 | h (moving platform vertical offset) | mm |
| 9 | L_tool (tool length) | mm |
| 10 | U-joint tilt angle | degrees |
| 11–15 | Base limb directions θ_b (5 values) | degrees |
| 16–20 | Moving limb directions θ_m (5 values) | degrees |
| 21–25 | Initial limb lengths l0 (5 values) | mm |

Parameters are loaded at runtime via:
```matlab
basic_paras = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'm');
```
Specifying `'unit', 'm'` converts mm (from Excel) to meters internally (`unit_para = 0.001`). Use `'unit', 'mm'` to keep mm.

## Key Conventions

### Variable Naming
- `p_seq` — 6×34 screw parameter matrix (exponential coordinates of initial transformations)
- `xi_seq` — 6×34 zero-configuration global screw coordinates
- `joint_q` — 6×5 joint variable matrix (one column per limb; SPR limb has 5 DOF, UPS limbs have 6)
- `T_ref` / `T_cal` — 4×4 homogeneous transformation matrices
- `Pos_ref_seq` — 5×N pose vector `[x; y; z; φ; θ]` (translations in m, angles in degrees)
- `B` — 3×5 base joint positions
- `P_m` — 3×5 moving platform joint positions (in platform frame)
- `limb_dir` — 5×2 matrix of limb direction angles `[θ_base, θ_move]` in radians

### Unit System
- Internal calculations: **meters, radians**
- Excel storage: **mm, degrees**
- Conversion handled by `basic_read()` via the `unit_para` factor

### Code Patterns
- All entry-point scripts begin with `clear` to reset workspace
- Iterative solvers use `err_max` tolerance (typically `1e-5` to `1e-9`) and `loop_max` safety limits
- Levenberg-Marquardt damping parameter `lambda` adapts based on gain ratio `eta` (calibration4.m)
- Total Least Squares via SVD of augmented matrix `[Jp_bar, -err]` (calibration5.m, `solve_tls.m`)
- Plots use Chinese labels with font `'微软雅黑'` and `'Times New Roman'`, figure size in centimeters
- Random seeds are set explicitly (e.g., `rng(0313+im)`) for reproducibility

### Screw Theory Convention
- Twist coordinates: `[ω; v]` (6×1), angular velocity followed by linear velocity
- Wrench coordinates: `[F; M]` (6×1), force followed by moment
- Reciprocal product operator: `Omega = [zeros(3,3) eye(3); eye(3) zeros(3,3)]`
- The SPR limb (limb 1) has joints: R-R-R-P-R (5 DOF)
- Each UPS limb (limbs 2–5) has joints: R-R-P-R-R-R (6 DOF)
- Joint axes: `zeta_r = [0;0;1;0;0;0]` (rotation about z), `zeta_p = [0;0;0;0;0;1]` (translation along z)

## Key Algorithms

### Forward Kinematics
Uses Newton-Raphson iteration on the body Jacobian. The closure constraint is that all five limbs must reach the same platform pose. The Jacobian `J_all` (24×24) relates passive joint velocities across limbs, and the update solves `J_all * Δq_passive = err`.

### Inverse Kinematics
Closed-form using Paden-Kahan subproblems:
- **PK1**: rotation about a single axis
- **PK2**: rotation about two intersecting axes
Active prismatic joints are solved directly from limb length; passive revolute joints use PK decomposition.

### Kinematic Calibration
Two methods are implemented:

**LM (calibration4.m, Chapter 3):**
1. Parameterize geometry into screw parameters `p_seq`
2. Compute row-space decomposition of the identification Jacobian
3. Iteratively solve `Jp_bar * Δp = err` with Levenberg-Marquardt damping
4. Rebuild `xi_seq` from updated `p_seq` each iteration
5. Recompute forward kinematics and pose error until convergence

**TLS (calibration5.m, Chapter 4):**
1. Same row-space decomposition framework (U, V_prep) to eliminate 92 redundant parameters
2. Form augmented matrix `C = [Jp_bar, -err]` (m × 113)
3. SVD of C; extract TLS solution `δp = -v(1:n) / v(n+1)` from last right singular vector
4. TLS simultaneously minimizes `||[ΔJ, Δr]||_F`, treating both Jacobian and residual as noisy
5. Fallback to pinv-LS when `σ_n ≈ σ_{n+1}` (degenerate case)
6. No damping parameter — regularization from SVD truncation

### Performance Indices
- **LTI** (Local Transmission Index) — local force/motion transmission quality
- **GCI** (Global Constraint Index) — constraint performance
- **OTI** (Output Transmission Index) — output motion transmission efficiency
- **ITI** (Input Transmission Index) — input transmission efficiency
These are computed by solving for TWS (Transmission Wrench Screw), OTS (Output Twist Screw), and SC (Constraint Screw) at each pose.

## Dependencies

- **MATLAB** (R2019b or later recommended — uses `readtable`, `deg2rad`, `surrogateopt`)
  - Optimization Toolbox (for `surrogateopt` in `parameter_optimize_surrogate.m`)
  - Statistics and Machine Learning Toolbox (used in some scripts)
  - No other toolboxes required for core kinematics
- **parameters.xlsx** must be present in the project root
- No external MATLAB toolboxes or MEX files
- The `lib_calib/` directory contains the TLS solver (`solve_tls.m`)

## Git Workflow

- Branch: `main`
- No CI/CD configured
- Commit messages are short English descriptions (e.g., "calib_start", "joint_angle_search", "para_opt_o1")
- `.gitignore` excludes `other/` and `sim_para_process/` directories (note: these are already tracked but future additions there are ignored)

## External Simulation Data

`sim_para_process/` contains scripts for processing RecurDyn multi-body dynamics simulation output (CSV files with/without gravity). This is used to validate analytical models against numerical simulation.
