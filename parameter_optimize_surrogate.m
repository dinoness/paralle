clear
path_add()

fprintf('>>>>> SPR-4UPS Structural Parameter Optimization >>>>>\n');
fprintf('>>>= start (%s) =<<<\n', string(datetime('now', 'Format', 'HH:mm:ss')));
fprintf('Using surrogateopt (Global Optimization Toolbox)\n\n');

%% 1. 优化设置
% 决策变量: [H_mm, h_mm, r1_mm, r2_mm, a_deg, b_deg]
lb = [-20, 60, 100,  60,  40,  40];
ub = [  5, 150, 150,  70,  50,  80];

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
fprintf('Bounds: H=[0,20], h=[10,50], r1=[90,130], r2=[60,80], a=[25,50]\n');
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
