# AGENTS.md — SPR-4UPS 并联机构分析与优化项目

> 本文件面向 AI 编程助手。若你从未接触过本项目，请先阅读此文件。

---

## 1. 项目概述

本项目是一个基于 MATLAB 的机器人学研究代码库，研究对象为 **SPR-4UPS 并联机构**（5 自由度：3 平移 + 2 转动）。

机构构型：
- 1 条 SPR 支链（球铰-移动副-转动副），提供 1 个约束
- 4 条 UPS 支链（万向铰-移动副-球铰），各提供 1 个主动驱动

主要研究内容：
- 基于螺旋理论（Screw Theory）与 Lie 群 / Lie 代数的运动学建模
- 运动学正解 / 逆解（指数积公式，POE）
- 工作空间离散搜索与最大内接圆柱提取
- 传递性能指标计算（ITI / OTI / LTI / GCI / LCI）
- 结构参数优化（ surrogateopt / 自定义搜索）
- 几何误差标定（Levenberg-Marquardt 迭代）
- 静力学分析与仿真验证

---

## 2. 技术栈与运行环境

- **语言**：MATLAB（R2020b 或更高版本推荐）
- **必需工具箱**：
  - Parallel Computing Toolbox（大量脚本使用 `parfor` 加速网格搜索）
  - Optimization Toolbox（`surrogateopt`, `fmincon` 等）
- **无外部包管理器**：不存在 `pyproject.toml`、`package.json` 或 `requirements.txt`。所有依赖为 MATLAB 内置函数或本项目自研函数。
- **操作系统**：Windows（开发环境），但代码本身跨平台。

---

## 3. 目录结构与代码组织

```
.
├── lib/                    # 核心运动学与标定算法库
│   ├── parameterize.m              # 由几何参数生成 6×34 旋量参数 p_seq 与零位全局旋量 xi_seq
│   ├── keni_sol_inverse.m          # 运动学逆解
│   ├── keni_sol_forward.m          # 运动学正解（Newton-Raphson / LM 迭代）
│   ├── keni_sol_forward_once.m     # 单次正解（给定关节量求位姿）
│   ├── jacobian_body.m             # 体坐标雅可比
│   ├── jacobian_space.m            # 空间坐标雅可比
│   ├── spr4ups_build_limb_twists.m # 构建各支链关节旋量
│   ├── build_line_wrench.m         # 构造线力旋量
│   ├── svd_nullspace.m             # 基于 SVD 的零空间计算
│   ├── pos2trans.m / trans2pos.m   # 位姿向量 <=> 4×4 齐次变换矩阵
│   ├── calib_iter_*.m              # 标定迭代相关的矩阵构造与恢复函数
│   └── ...
├── lib_math/               # 螺旋理论与 Lie 代数基础数学库
│   ├── exp_se3.m / log_se3.m       # SE(3) 指数映射与对数映射
│   ├── exp_so3.m / log_so3.m       # SO(3) 指数映射与对数映射
│   ├── adjoint_m.m                 # 伴随矩阵
│   ├── screw_efficiency.m          # 螺旋效率（归一化互易积）
│   ├── screw_apparent_power.m      # 螺旋视在功率
│   ├── screw_unitize.m             # 螺旋单位化
│   ├── recip.m                     # 互易积
│   ├── Paden_Kahan1/2.m            # Paden-Kahan 子问题求解
│   └── skew.m                      # 反对称矩阵
├── lib_para/               # 参数读取
│   └── basic_read.m                # 从 parameters.xlsx 读取几何参数，支持列选择与单位转换
├── lib_opt_struct/         # 结构优化与性能指标
│   ├── evaluate_spr4ups_objective.m # surrogateopt 目标函数（OTI 均值 + 5% 分位数）
│   └── compute_ltigci.m             # 计算 LTI / GCI / OTI / ITI
├── lib_calib/              # 标定专用库（当前为空或待扩展）
├── AI_code/                # AI 辅助生成的代码片段
│   └── SPR_4UPS_Statics.m          # 基于力雅可比的静力学求解（给定外力求驱动力与约束力）
├── other/                  # 历史脚本、可视化辅助、测试动画
│   ├── ArcDraw.m, Circle.m         # 几何绘图类
│   ├── Workspace.m, Workspace2.m   # 旧版工作空间脚本
│   ├── kinematics_test.m           # 运动学测试
│   └── *.gif                       # 生成的演示动画
├── sim_para_process/       # 与 RecurDyn 仿真对比数据
│   ├── recurdyn0427.m
│   ├── no_gravity.csv
│   └── with_gravity.csv
├── parameters.xlsx         # 中央参数文件（多列对应不同工况/优化阶段）
├── optimization_result.mat # 结构优化结果存档
└── [根目录脚本]            # 见下文“主要入口脚本”
```

### 主要入口脚本（根目录）

| 脚本名 | 用途 |
|--------|------|
| `workspace_discrete_v3.m` | **核心工作空间分析**：离散网格搜索、含姿态可达空间、最大内接圆柱提取、LTI/GCI/LCI 计算、可视化。大量使用 `parfor`。 |
| `Joint_angle_search.m` | 球铰链最优倾斜角 `phi` 搜索（以圆柱体积为目标），分多阶段粗搜+细搜。 |
| `parameter_optimize.m` / `parameter_optimize2.m` / `parameter_optimize_surrogate.m` | 结构参数优化入口，分别对应不同优化策略。 |
| `calibration2.m` ~ `calibration4.m` | 几何误差标定脚本，迭代辨识结构参数。`calibration4.m` 为最新版，使用 LM 阻尼与行空间分解。 |
| `velocity_static_force.m` | 静力学与速度映射分析，绘制力矢量图。 |
| `force_analysis_screw.m` | 基于螺旋理论的 ITI/OTI/LTI/GCI 计算与作图。 |
| `structure_err.m` | 给定结构误差（如基座铰链点偏移），观察位姿变化。 |
| `keni_sol.m` | 运动学正逆解的演示与验证脚本。 |
| `scerw_theory.m` / `screw_verify.m` | 螺旋理论验证与自由度分析。 |
| `test.m` / `test_script.m` | 临时测试与 AI 生成代码调用示例。 |
| `path_add.m` | **环境初始化**：将 `lib`, `lib_para`, `lib_math`, `lib_opt_struct`, `lib_calib` 加入 MATLAB 路径。 |

---

## 4. 关键开发约定

### 4.1 单位约定（极易出错）

- 全局变量 `unit_para` 控制单位换算：
  - `unit_para = 0.001` 表示程序内部使用 **米 (m)**，而 Excel 中参数以 **毫米 (mm)** 记录。
  - `unit_para = 1` 表示程序内部使用 **毫米 (mm)**。
- `basic_read.m` 支持 `'unit', 'm'` 或 `'unit', 'mm'` 参数，自动完成转换。
- **角度**：大部分脚本中，位姿向量的角度单位是 **度 (deg)**，但在调用 `pos2trans` 时需显式指定 `'unit', 'rad'`  if 输入为弧度。

### 4.2 参数文件 `parameters.xlsx`

- 单文件多列管理不同配置：
  - 列 `B`：默认 / 基础参数
  - 列 `C`、`D`、`M` 等：不同优化阶段、标定实验或工作空间分析用的参数集。
- 读取方式：`basic_read('parameters.xlsx', 'column', 'M', 'unit', 'm')`。
- 修改参数前务必确认当前脚本使用的是哪一列。

### 4.3 MATLAB 路径管理

- 每个主脚本开头通常调用 `path_add()` 或 `addpath(genpath('./lib'))`。
- **不要在脚本中永久保存路径**（避免 `savepath`），保持项目自包含。

### 4.4 并行计算

- 工作空间搜索大量使用 `parfor`。运行前建议：
  ```matlab
  parpool('local', 5);  % 根据 CPU 核心数调整
  ```
- `parfor` 循环内部变量必须满足切片变量规则；合并结果通常使用 `cell` + 事后拼接。

### 4.5 代码风格

- **注释语言**：中文为主，部分函数头部使用英文文档字符串。
- **变量命名**：
  - 蛇形命名与驼峰命名混用，整体偏向小写+下划线。
  - 常见缩写：`seq`（序列）、`limb`（支链）、`plant`（平台）、`pos`（位置）、`T_ref`（参考位姿变换矩阵）、`p_seq`（旋量参数矩阵）。
- **图形输出**：
  - 统一白色背景：`figure('Color', [1 1 1])`
  - 常用字体：`Times New Roman`（坐标轴）与 `微软雅黑`（标签）。
  - 固定配色：`#FF7F50`（基座/珊瑚色）、`#32CD32`（动平台/绿色）、`#4682B4`（支链/钢蓝色）。

### 4.6 版本控制习惯

- 存在大量同名不同版本的文件（如 `workspace_discrete.m`, `v2.m`, `v3.m`）。
- 旧版本通常保留在根目录或 `other/` 中，**修改前请确认正在编辑的是最新版本**（通常是版本号最大的）。
- `.gitignore` 忽略了 `/other` 和 `/sim_para_process`，这两个目录下的内容不会被 Git 跟踪。

---

## 5. 运行与测试

### 5.1 如何运行一个典型分析

以工作空间分析为例：

```matlab
% 1. 进入项目根目录
% 2. 在 MATLAB 命令行执行：
path_add();
workspace_discrete_v3
```

脚本会自动：
1. 从 `parameters.xlsx` 读取参数；
2. 执行离散网格搜索（`parfor` 加速）；
3. 提取最大内接圆柱；
4. 计算 LTI / GCI / LCI；
5. 绘制 3D 散点图与机构简图。

### 5.2 如何运行结构优化

```matlab
parameter_optimize_surrogate   % 使用 surrogateopt 进行全局优化
```

或：

```matlab
% 手动调用目标函数
J = evaluate_spr4ups_objective([H_mm, h_mm, r1_mm, r2_mm, a_deg, b_deg], 1.0, 1.0);
```

### 5.3 如何验证运动学

```matlab
keni_sol          % 包含正逆解演示
screw_verify      % 螺旋理论验证
```

### 5.4 测试策略

- **无单元测试框架**：本项目为科研代码，未使用 MATLAB Unit Test Framework。
- 验证方式：
  - 运行 `screw_verify.m`、`verify_OTI_3prs.m`、`verify_OTI_6ups.m` 等对比脚本；
  - 将 MATLAB 计算结果与 `sim_para_process/recurdyn0427.m` 中的仿真数据对比；
  - 检查 `parameterize.m` 生成的 `p_seq` 与手动推导是否一致。
- 修改核心库（`lib/` 或 `lib_math/`）后，建议依次运行：
  1. `keni_sol.m` —— 正逆解闭环验证
  2. `screw_verify.m` —— 螺旋理论一致性
  3. `workspace_discrete_v3.m` —— 端到端工作空间分析（取极小网格快速验证）

---

## 6. 安全与数值稳定性注意事项

- **奇异位形**：在 Home Position 附近，SPR 支链的 R 关节轴线可能与支链方向平行，导致力雅可比矩阵 `G_matrix` 奇异。`SPR_4UPS_Statics.m` 与 `velocity_static_force.m` 中已加入 `cond` 检查与 `warning`。
- **标定迭代发散**：`calibration4.m` 使用自适应阻尼 LM 方法，但初始值离真值过远时仍可能不收敛。建议先用小扰动测试。
- **parfor 中的大矩阵拼接**：`workspace_discrete_v3.m` 等脚本在 `parfor` 中使用动态数组拼接（`[work_space_ang local_ws_ang]`），这在 MATLAB 中可行但通信开销较大。修改时尽量避免在 `parfor` 内频繁改变大矩阵尺寸。
- **持久变量（persistent）**：`evaluate_spr4ups_objective.m` 使用 `persistent basic_paras_base` 缓存参数读取结果。若 `parameters.xlsx` 内容在 MATLAB 会话期间被修改，需执行 `clear evaluate_spr4ups_objective` 才能重新加载。

---

## 7. 快速参考：常用函数接口

```matlab
% 参数化
[p_seq, xi_seq] = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

% 逆解
joint_q = keni_sol_inverse(T_ref, B, l0_seq, P_m, p_seq);

% 正解（迭代）
[T, joint_q] = keni_sol_forward(joint_q0, p_seq, err_max);

% 位姿向量 <-> 变换矩阵
T = pos2trans([x; y; z; phi; theta], B, 'unit', 'rad');
pos = trans2pos(T);

% 性能指标
[LTI, GCI, OTI, ITI] = compute_ltigci(T_ref, B, l0_seq, P_m, p_seq);

% 参数读取
paras = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'm');
```

---

## 8. 扩展与修改建议

- **新增标定算法**：建议放入 `lib_calib/`，并在根目录新建 `calibrationX.m` 作为入口。
- **新增性能指标**：在 `lib_opt_struct/` 中新增函数，保持输入接口与 `compute_ltigci.m` 一致（`T_ref, B, l0_seq, P_m, p_seq`）。
- **新增数学工具**：在 `lib_math/` 中新增，确保函数头部注释包含输入输出维度说明。
- **修改 `parameters.xlsx` 列结构时**：务必同步更新 `lib_para/basic_read.m` 的读取范围（`range_start`, `range_end`）与 `table2array` 列映射逻辑。

---

*最后更新：2026-05-29*
