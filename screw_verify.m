% 先运行force_analysis_screw.m

%% ========== 并联机构旋量验证（修正版） ==========
CONVENTION = 'A';  % Twist=[w;v], Wrench=[f;m]
Delta = [zeros(3), eye(3); eye(3), zeros(3)];

fprintf('=================================================\n');
fprintf('并联机构旋量计算结果验证（修正版）\n');
fprintf('=================================================\n');

%% 1. SI 输入关节旋量验证（P副检查）
fprintf('\n【1. SI 输入关节旋量验证】\n');
fprintf('SI 是 P 副旋量，前3维（角速度）应≈0，rank 应≤3\n');
for i = 1:5
    si = SI(:,i);
    omega = si(1:3); v = si(4:6);
    fprintf('SI(%d): |ω|=%.2e, |v|=%.4f, v_dir=[%5.2f,%5.2f,%5.2f]\n', ...
        i, norm(omega), norm(v), v(1), v(2), v(3));
    assert(norm(omega) < 1e-6, '❌ SI(%d) 不是纯平移！P副旋量前3维必须为0', i);
end
r_si = rank(SI);
fprintf('rank(SI) = %d (期望: ≤3)\n', r_si);
assert(r_si <= 3, '❌ P副旋量矩阵秩异常（>3说明不是纯平移）');

%% 2. SO 输出运动旋量验证
fprintf('\n【2. SO 输出运动旋量验证】\n');
r_so = rank(SO);
fprintf('rank(SO) = %d (期望: 5)\n', r_so);
assert(r_so == 5, '❌ SO 不满秩，机构处于奇异位形或计算错误');

%% 3. 互易关系验证（核心：教材公式 9.1-13）
fprintf('\n【3. 互易关系验证】\n');

% 3.1 约束力 SC 与所有输出运动 SO 互易
recip_SC_SO = SC' * Delta * SO;
fprintf('SC·SO (应全≈0): ');
fprintf('%.2e  ', recip_SC_SO); fprintf('\n');
assert(all(abs(recip_SC_SO) < 1e-8), '❌ SC 与 SO 不互易');

% 3.2 传递力 ST(j) 与输出运动 SO(i) 互易 (j≠i) —— 教材公式 9.1-13
fprintf('\nST(j)·SO(i) (j≠i, 应≈0):\n');
for i = 1:5
    for j = 1:5
        if j ~= i
            val = ST(:,j)' * Delta * SO(:,i);
            fprintf('  ST(%d)·SO(%d) = %.2e  ', j, i, val);
            if abs(val) > 1e-6
                fprintf('❌\n');
                error('ST(%d) 与 SO(%d) 应互易 (j≠i)', j, i);
            else
                fprintf('✓\n');
            end
        end
    end
end

% 3.3 同分支力与运动可做功 (j=i)
fprintf('\nST(i)·SO(i) (应≠0，表示能做功):\n');
for i = 1:5
    val = ST(:,i)' * Delta * SO(:,i);
    fprintf('  ST(%d)·SO(%d) = %.4f  ', i, i, val);
    if abs(val) < 1e-6
        fprintf('❌\n');
        warning('ST(%d) 与 SO(%d) 互易积为0，该驱动无法做功', i, i);
    else
        fprintf('✓\n');
    end
end

%% 4. 力空间完整性
fprintf('\n【4. 力空间完整性】\n');
r_f = rank([ST, SC]);
fprintf('rank([ST,SC]) = %d (期望: 6)\n', r_f);
assert(r_f == 6, '❌ 力旋量空间不满秩（机构奇异或 ST/SC 计算错误）');

%% 5. SI 与 SO 的映射关系（如有分支雅可比）
% 如果你保存了各分支的被动关节旋量矩阵，可以验证：
% SO(:,i) 应该落在第 i 个分支的关节旋量空间内（包含主动+被动关节）
fprintf('\n【5. 物理合理性检查】\n');
fprintf('若你记录了分支几何参数，请手动检查：\n');
fprintf('  • UPS 分支的 ST(:,i) 轴线是否通过对应球副中心 S_i\n');
fprintf('  • SPR 分支的 ST(:,5) 是否沿其 P 副轴线方向\n');
fprintf('  • SO 的各列是否满足 SC 的约束（已自动验证）\n');

fprintf('\n=================================================\n');
fprintf('所有验证通过！\n');
fprintf('=================================================\n');