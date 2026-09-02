# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MATLAB codebase for kinematic modeling, calibration, and workspace analysis of an **SPR-4UPS parallel mechanism** (5-limb: 1 SPR + 4 UPS). Uses **screw theory** and the **Product of Exponentials (POE)** formulation throughout. Reference paper: `c02 并联机构的运动学误差建模及参数可辨识性分析_孔令雨.pdf`.

## Essential Setup

Every script must call `path_add()` before using library functions:

```matlab
path_add();  % adds lib/, lib_para/, lib_math/, lib_opt_struct/, lib_calib/ to path
```

All geometric parameters are read from `parameters.xlsx` via `basic_read()`:

```matlab
basic_paras = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'm');
% column: which Excel column to read (B, C, M, etc.)
% unit: 'm' → values in meters (multiplies by 0.001); 'mm' → values in mm (raw)
```

## Library Architecture

| Directory | Purpose |
|---|---|
| `lib/` | Core kinematics: POE parameterization, forward/inverse kinematics, body/space Jacobians, calibration matrices |
| `lib_math/` | Lie group operations: `exp_se3`/`log_se3`, `exp_so3`/`log_so3`, Paden-Kahan subproblems, screw reciprocity/power/efficiency |
| `lib_calib/` | Calibration pipeline: measurement sequence generation (`calib_seq_generate`), TLS solver (`solve_tls`), 3-point transform (`three_pts2trans`) |
| `lib_opt_struct/` | Structural parameter optimization: `compute_ltigci` (OTI metric), `evaluate_spr4ups_objective` |
| `lib_para/` | `basic_read` — reads nominal geometric parameters from Excel |

## Key Data Structures

- **`p_seq`** (6×34 matrix): Kinematic parameters in screw coordinates. Columns 1–6 = SPR limb, 7–13/14–20/21–27/28–34 = four UPS limbs. Each column is the exponential coordinate (se(3) twist) of a joint-to-joint transform. This is the central parameter representation — calibration iterates on `p_seq`.
- **`xi_seq`** (6×34 matrix): Zero-configuration joint twist coordinates in the global frame. Rebuilt from `p_seq` after each calibration iteration via `rebuild_xi_seq_from_p`.
- **`joint_q`** (6×5 matrix): Joint displacements. Column = limb (1–5), row = joint axis within that limb. SPR limb has an extra passive joint.
- **`Pos_ref_seq`** (5×n matrix): Reference poses as `[x; y; z; theta_y; theta_z]` in SI units (m, rad).
- **`basic_paras`** struct: Nominal geometry including `R1, R2` (base radii), `r1, r2` (platform radii), `H, h` (offsets), `B` (base joint positions), `P_m` (platform joint positions), `l0_seq` (initial limb lengths), `limb_dir` (joint direction angles).

## Parameterization Pipeline

```
Geometric parameters (B, P_m, r1, r2, limb_dir, l0_seq, joint_u_angle_tilt)
    │
    ▼ parameterize()          ──── produces p_seq (6×34), xi_seq (6×34)
    │
    ▼ keni_sol_inverse()      ──── joint_q from desired pose
    │
    ▼ keni_sol_forward()      ──── actual pose from joint_q (Newton iterative)
    │
    ▼ deparameterize()        ──── extract geometric params back from p_seq
```

`parameterize()` and `deparameterize()` are inverses. `parameterize` encodes each joint-to-joint transformation as an se(3) twist in `p_seq`. For the SPR limb, all joint axes are derived analytically; for UPS limbs, the U-joint axis direction depends on `joint_u_angle_tilt` and `limb_dir`.

## Forward Kinematics (keni_sol_forward)

Iterative Newton-type solver. Given joint displacements and current `p_seq`:
1. `keni_sol_forward_once()` computes an initial pose guess by composing exponentials along all limbs, averaging the limb chain results.
2. Body Jacobians are computed, active joints are stripped, and the coupled system `J_all * Δjoint_passive = err` is solved.
3. Passive joint values are updated, pose is recomputed, and the loop repeats until the pose error across limbs drops below `err_max`.

## Calibration Workflow

The calibration scripts (`calibration3.m`, `calibration4.m`, `calibration5.m`) follow the pattern in Chapter 3/4 of the reference paper:

1. **Generate nominal model**: `parameterize()` from nominal geometry
2. **Create "true" model**: Apply systematic bias (`delta_B_sys`) + random parametric errors to simulate manufacturing imperfections
3. **Generate measurements**: Nominal inverse kinematics → joint commands → "true" forward kinematics → add Gaussian measurement noise
4. **Iterative identification**: Levenberg-Marquardt (LM) in `calibration4.m` or Regularized TLS in `calibration5.m`
   - `calib_iter_row_space_matrix()` — computes nullspace/redundancy structure
   - `calib_iter_matrix2()` — builds the identification Jacobian `Jp_bar`
   - `restore_full_param_increment()` — maps parameter increments back through the nullspace structure
   - Adaptive damping factor λ updated by gain ratio η
5. **Validation**: Compare position/orientation errors on a holdout validation set before and after calibration
6. **Export**: `deparameterize()` → write `calibrated_params.csv`

## Main Script Reference

| Script | Purpose |
|---|---|
| `calibration4.m` | LM-based kinematic calibration (primary, most current) |
| `calibration5.m` | Regularized TLS calibration (alternative to LM, Chapter 4) |
| `calibration3.m` | Earlier POE calibration experiment (deprecated) |
| `workspace_discrete_v3.m` | Workspace search with OTI/LCI computation, cylinder space analysis |
| `Joint_angle_search.m` | Joint angle limit search across the workspace |
| `force_analysis_screw.m` | Screw-theory based static force analysis |
| `structure_err.m` | Sensitivity analysis — how structural parameter errors affect pose |
| `parameter_optimize_surrogate.m` | Structural parameter optimization via `surrogateopt` |
| `keni_sol.m` | Standalone forward kinematics demonstration |
| `test_script.m` | Ad-hoc testing scratchpad |
| `controller_verify.m` | Controller verification routine |
| `screw_verify.m` | Screw theory verification |
| `velocity_static_force.m` | Velocity and static force analysis |
| `test_deparameterize.m` | Round-trip test of parameterize → deparameterize |

## C Interface

`read_calibrated_params.h` is a standalone C header that reads `calibrated_params.csv` into a `CalibratedParams` struct. It's a portable way to consume calibrated parameters in C-based controllers. The CSV uses `#` for comments and label-prefix conventions (`r1`, `r2`, `alpha`, `l0_N`, `B_N`, `Pm_N`, `limb_dir_N`). All values are SI (m, rad).

## Data Files

- `parameters.xlsx` — Nominal geometric parameters. Columns B/C/M/etc. hold different parameter sets. Row layout: l_max, l_min, R1, R2, H, r1, r2, h, L_tool, joint_u_angle_tilt, then 5× limb_dir base angles, 5× limb_dir move angles, 5× l0 values.
- `calibrated_params.csv` — Calibrated output (SI units: m, rad). 12 decimal places. `#`-prefixed header lines.
- `optimization_result.mat` — Saved `surrogateopt` results from structural optimization runs.
- `sim_para_process/` — RecurDyn simulation data (CSV) for verifying gravity/no-gravity cases.
