clear
path_add()

fprintf('>>>>> SPR-4UPS Structural Parameter Optimization >>>>>\n');
fprintf('>>>= start (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));
fprintf('Using surrogateopt (Global Optimization Toolbox)\n\n');

%% 1. 优化设置
% 决策变量: [H_mm, h_mm, r1_mm, r2_mm, a_deg, b_deg]
lb = [-20,  80, 100,  60,  40,  50];
ub = [ 30, 200, 180,  70,  50,  70];

% 权重（可调整）
k1 = 0.6;
k2 = 0.4;

% 目标函数句柄
fun = @(x) evaluate_spr4ups_objective(x, k1, k2);

% surrogateopt 选项
options = optimoptions('surrogateopt', ...
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 500, ...
    'PlotFcn', 'surrogateoptplot', ...
    'UseParallel', false);   % parfor 已在目标函数内部使用

rng(42);  % 可重复性

%% 2. 运行优化
fprintf('Starting surrogateopt...\n');
fprintf('Bounds: H=[%.0f,%.0f], h=[%.0f,%.0f], r1=[%.0f,%.0f], r2=[%.0f,%.0f], a=[%.0f,%.0f], b=[%.0f,%.0f]\n', ...
    lb(1), ub(1), lb(2), ub(2), lb(3), ub(3), lb(4), ub(4), lb(5), ub(5), lb(6), ub(6));
fprintf('Weights: k1=%.2f, k2=%.2f\n\n', k1, k2);

[x_opt, fval, exitflag, output] = surrogateopt(fun, lb, ub, options);

%% 3. 输出优化结果
fprintf('\n========== Optimization Result (Coarse Grid) ==========\n');
fprintf('H  = %.4f mm\n', x_opt(1));
fprintf('h  = %.4f mm\n', x_opt(2));
fprintf('r1 = %.4f mm\n', x_opt(3));
fprintf('r2 = %.4f mm\n', x_opt(4));
fprintf('a  = %.4f deg\n', x_opt(5));
fprintf('b  = %.4f deg\n', x_opt(6));
fprintf('J  = %.6f\n', fval);
fprintf('Function evaluations: %d\n', output.funccount);
fprintf('=======================================================\n\n');

%% 4. 保存结果
results.x_opt = x_opt;
results.fval = fval;
results.k1 = k1;
results.k2 = k2;
results.lb = lb;
results.ub = ub;
results.output = output;
save('optimization_result.mat', 'results');
fprintf('Results saved to optimization_result.mat\n');

fprintf('>>>= done (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));
