% test_deparameterize.m — 验证 deparameterize 正确性（parameterize 逆映射）
clear
path_add();
fprintf('===== deparameterize 正确性测试 =====\n');

tol = 1e-12;

%% 1. 加载名义参数
basic_paras = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'm');
B_nom      = basic_paras.B;
r1_nom     = basic_paras.r1;
r2_nom     = basic_paras.r2;
l0_nom     = basic_paras.l0_seq;
P_m_nom    = basic_paras.P_m;
limb_nom   = basic_paras.limb_dir;
alpha_nom  = basic_paras.joint_u_angle_tilt;

%% 2. parameterize → deparameterize 往返测试
[p_seq, ~, info] = parameterize(limb_nom, B_nom, r1_nom, r2_nom, l0_nom, P_m_nom, alpha_nom);
[B_out, r1_out, r2_out, l0_out, P_m_out, limb_out, alpha_out] = ...
    deparameterize(p_seq, P_m_nom, limb_nom);

%% 3. 逐项比对
n_pass = 0;
n_fail = 0;

fprintf('\n--- 基座铰点 B (3×5) ---\n');
for i = 1:5
    err = norm(B_out(:,i) - B_nom(:,i));
    if err < tol, n_pass = n_pass + 1;
    else, fprintf('  ✗ B(:,%d) 误差 = %.2e\n', i, err); n_fail = n_fail + 1;
    end
end
if ~any(n_fail); fprintf('  ✓ 全部通过\n'); end

fprintf('\n--- 平台半径 ---\n');
if abs(r1_out - r1_nom) < tol
    fprintf('  ✓ r1 = %.6f (误差 %.2e)\n', r1_out, abs(r1_out - r1_nom)); n_pass = n_pass + 1;
else
    fprintf('  ✗ r1 = %.6f (期望 %.6f)\n', r1_out, r1_nom); n_fail = n_fail + 1;
end
if abs(r2_out - r2_nom) < tol
    fprintf('  ✓ r2 = %.6f (误差 %.2e)\n', r2_out, abs(r2_out - r2_nom)); n_pass = n_pass + 1;
else
    fprintf('  ✗ r2 = %.6f (期望 %.6f)\n', r2_out, r2_nom); n_fail = n_fail + 1;
end

fprintf('\n--- 初始杆长 l0_seq (1×5) ---\n');
for i = 1:5
    err = abs(l0_out(i) - l0_nom(i));
    if err < tol, n_pass = n_pass + 1;
    else, fprintf('  ✗ l0(%d) 误差 = %.2e\n', i, err); n_fail = n_fail + 1;
    end
end
if ~any(n_fail); fprintf('  ✓ 全部通过\n'); end

fprintf('\n--- 动平台铰点 P_m(3,:) ---\n');
for i = 1:5
    err = abs(P_m_out(3,i) - P_m_nom(3,i));
    if err < tol, n_pass = n_pass + 1;
    else, fprintf('  ✗ P_m(3,%d) 误差 = %.2e\n', i, err); n_fail = n_fail + 1;
    end
end
if ~any(n_fail); fprintf('  ✓ 全部通过\n'); end

fprintf('\n--- 支链方向角 limb_dir(2:5,:) ---\n');
wrap = @(a) mod(a + pi, 2*pi) - pi;  % 包裹到 [-π, π]
for i = 2:5
    d1 = wrap(limb_out(i,1) - limb_nom(i,1));
    d2 = wrap(limb_out(i,2) - limb_nom(i,2));
    if abs(d1) < tol, n_pass = n_pass + 1;
    else, fprintf('  ✗ limb(%d,1) 误差 = %.2e rad\n', i, d1); n_fail = n_fail + 1;
    end
    if abs(d2) < tol, n_pass = n_pass + 1;
    else, fprintf('  ✗ limb(%d,2) 误差 = %.2e rad\n', i, d2); n_fail = n_fail + 1;
    end
end
if ~any(n_fail); fprintf('  ✓ 全部通过\n'); end

fprintf('\n--- U 副倾角 ---\n');
if abs(alpha_out - alpha_nom) < tol
    fprintf('  ✓ alpha = %.6f rad (%.4f°) (误差 %.2e)\n', alpha_out, rad2deg(alpha_out), abs(alpha_out - alpha_nom));
    n_pass = n_pass + 1;
else
    fprintf('  ✗ alpha 误差 = %.2e rad\n', abs(alpha_out - alpha_nom)); n_fail = n_fail + 1;
end

%% 4. 往返一致性：再 parameterize 一次看 p_seq 是否一致
fprintf('\n--- 往返闭合性 (p_seq 再生成) ---\n');
[p_seq2, ~] = parameterize(limb_out, B_out, r1_out, r2_out, l0_out, P_m_out, alpha_out);
p_err = max(abs(p_seq2(:) - p_seq(:)));
if p_err < tol
    fprintf('  ✓ p_seq 往返闭合, max|Δ| = %.2e\n', p_err);
    n_pass = n_pass + 1;
else
    fprintf('  ✗ p_seq 往返闭合失败, max|Δ| = %.2e\n', p_err); n_fail = n_fail + 1;
end

%% 5. 结果汇总
fprintf('\n===== 结果: %d 通过, %d 失败 =====\n', n_pass, n_fail);
if n_fail == 0
    fprintf('deparameterize 正确性验证通过。\n');
end
