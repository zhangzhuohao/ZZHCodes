%% ======================== 完整代码：无外部时间信号的双模块RNN模型 ========================
clear; close all; clc;

%% 1. 参数设置
% 网络结构
n_m2 = 50;              % M2模块神经元数
n_down = 20;            % 下游模块神经元数
inputSize = 2;          % 输入: [任务指令1, 任务指令2] (已移除时间信号!)
outputSize = 1;         % 输出: 运动命令

% 时间参数
totalLength = 150;      % 试次总长度 = baselineLength + taskLength
baselineLength = 30;    % 基线期长度
taskLength = totalLength - baselineLength; % 任务期长度

% 训练参数
numEpochs = 500;        % 增加训练轮数（内部计时需要更多训练）
batchSize = 32;         % 批大小
learningRate = 0.001;   % 学习率
cueProb = 0.5;          % 线索触发任务概率

fprintf('无外部时间信号模型参数:\n');
fprintf('  输入维度: %d (无外部时间信号)\n', inputSize);
fprintf('  网络结构: M2(%d) -> 下游(%d)\n', n_m2, n_down);
fprintf('  时间结构: 总长度=%d (基线%d + 任务%d)\n', totalLength, baselineLength, taskLength);
fprintf('  训练: %d轮, 批大小=%d\n', numEpochs, batchSize);

%% 3. 初始化网络参数 (Xavier初始化)
fprintf('\n初始化网络参数...\n');

% M2模块参数
W_m2_rec = randn(n_m2, n_m2) * sqrt(2 / (n_m2 + n_m2));
W_m2_in = randn(n_m2, inputSize) * sqrt(2 / (n_m2 + inputSize));  % 注意: inputSize现在是2
b_m2 = zeros(n_m2, 1);

% M2 -> 下游连接
W_m2_to_down = randn(n_down, n_m2) * sqrt(2 / (n_down + n_m2));

% 下游模块参数
W_down_rec = randn(n_down, n_down) * sqrt(2 / (n_down + n_down));
b_down = zeros(n_down, 1);

% 输出层参数
W_out = randn(outputSize, n_down) * sqrt(2 / (outputSize + n_down));
b_out = zeros(outputSize, 1);

% 打包参数
params = struct('W_m2_rec', W_m2_rec, 'W_m2_in', W_m2_in, 'b_m2', b_m2, ...
    'W_m2_to_down', W_m2_to_down, ...
    'W_down_rec', W_down_rec, 'b_down', b_down, ...
    'W_out', W_out, 'b_out', b_out);

%% 5. 训练循环
fprintf('\n开始训练...\n');
lossHistory = zeros(1, numEpochs);
bestLoss = inf;
bestParams = params;

for epoch = 1:numEpochs
    % 生成训练数据
    [inputs, targets, trialInfo] = generateBatch(batchSize, totalLength, baselineLength, cueProb);

    % 前向传播和反向传播
    [loss, grads, ~, ~, ~] = forwardBackward(params, inputs, targets, baselineLength);
    lossHistory(epoch) = loss;

    % 梯度裁剪
    clipValue = 1.0;
    for field = fieldnames(grads)'
        grads.(field{1}) = sign(grads.(field{1})) .* min(abs(grads.(field{1})), clipValue);
    end

    % 更新参数
    for field = fieldnames(params)'
        params.(field{1}) = params.(field{1}) - learningRate * grads.(field{1});
    end

    % 保存最佳参数
    if loss < bestLoss
        bestLoss = loss;
        bestParams = params;
    end

    % 每50轮显示进度
    if mod(epoch, 50) == 0
        fprintf('  轮次 %d/%d, 损失: %.6f\n', epoch, numEpochs, loss);
    end
end

fprintf('训练完成! 最佳损失: %.6f\n', bestLoss);
params = bestParams;

%% 6. 测试与可视化
fprintf('\n测试模型...\n');
testBatchSize = 10;

% 生成测试数据
[testInputs, testTargets, testInfo] = generateBatch(testBatchSize, totalLength, baselineLength, 0.5);

% 前向传播
[~, ~, h_m2_all_test, h_down_all_test, testOutputs] = forwardBackward(params, testInputs, testTargets, baselineLength);

% 计算测试损失
testLoss = 0;
taskStepsCount = 0;
for t = baselineLength+1:totalLength
    y_t = reshape(testOutputs(:, t, :), [], testBatchSize);
    target_t = reshape(testTargets(:, t, :), [], testBatchSize);
    testLoss = testLoss + sum((y_t(:) - target_t(:)).^2);
    taskStepsCount = taskStepsCount + 1;
end
testLoss = testLoss / (testBatchSize * taskStepsCount);
fprintf('测试损失: %.6f\n', testLoss);

%% 7. 可视化结果
figure('Position', [100, 100, 1400, 800]);

% 1. 训练损失曲线
subplot(3, 4, [1, 2]);
plot(lossHistory, 'LineWidth', 2);
xlabel('训练轮次'); ylabel('损失');
title('训练损失曲线 (无外部时间信号)');
grid on;
xlim([1, numEpochs]);

% 2. 示例试次：输入、目标和模型输出
subplot(3, 4, [3, 4]);
trialIdx = 1;
time = 1:totalLength;

% 绘制任务指令
plot(time, squeeze(testInputs(trialIdx, :, 1)), 'b-', 'LineWidth', 1, 'DisplayName', '任务指令1(线索)');
hold on;
plot(time, squeeze(testInputs(trialIdx, :, 2)), 'r-', 'LineWidth', 1, 'DisplayName', '任务指令2(自主)');
plot(time, squeeze(testTargets(trialIdx, :, 1)), 'k-', 'LineWidth', 2, 'DisplayName', '目标运动');
plot(time, squeeze(testOutputs(trialIdx, :, 1)), 'g--', 'LineWidth', 1.5, 'DisplayName', '模型输出');

% 标记基线期和任务期
ylim = get(gca, 'YLim');
fill([1, baselineLength, baselineLength, 1], [ylim(1), ylim(1), ylim(2), ylim(2)], ...
    [0.9, 0.9, 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
text(baselineLength/2, ylim(2)*0.9, '基线期', 'HorizontalAlignment', 'center');
text((baselineLength+totalLength)/2, ylim(2)*0.9, '任务期', 'HorizontalAlignment', 'center');

xlabel('时间步'); ylabel('强度');
title(sprintf('示例试次 (任务类型: %s)', testInfo(trialIdx).isCued));
legend('Location', 'best');
grid on;
hold off;

% 3. M2群体活动 (热图，前10个神经元)
subplot(3, 4, [5, 6]);
numNeuronsToShow = min(50, n_m2);
m2_activity = squeeze(mean(h_m2_all_test(:, :, 1:numNeuronsToShow), 1))';
imagesc(time, 1:numNeuronsToShow, m2_activity);
xlabel('时间步'); ylabel('神经元索引');
title('M2群体活动 (平均)');
colorbar;
colormap('jet');

hold on;
plot([baselineLength, baselineLength], [0.5, numNeuronsToShow+0.5], 'w--', 'LineWidth', 1.5);
hold off;

% 4. 任务类型对比：M2平均活动
subplot(3, 4, [7, 8]);

% 分离两种任务
cuedTrials = find([testInfo.isCued]);
selfTrials = find(~[testInfo.isCued]);

if ~isempty(cuedTrials)
    m2_mean_cued = mean(mean(h_m2_all_test(cuedTrials, :, :), 3), 1);
    plot(time, m2_mean_cued, 'b-', 'LineWidth', 2, 'DisplayName', '线索触发任务');
    hold on;
end

if ~isempty(selfTrials)
    m2_mean_self = mean(mean(h_m2_all_test(selfTrials, :, :), 3), 1);
    plot(time, m2_mean_self, 'r-', 'LineWidth', 2, 'DisplayName', '自主计时任务');
end

% 标记基线期
ylim = get(gca, 'YLim');
fill([1, baselineLength, baselineLength, 1], [ylim(1), ylim(1), ylim(2), ylim(2)], ...
    [0.9, 0.9, 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlabel('时间步'); ylabel('平均活动');
title('不同任务下M2平均活动对比 (无外部时间信号)');
legend('Location', 'best');
grid on;
hold off;

% 5. 状态空间轨迹 (PCA降维)
subplot(3, 4, [9, 10]);

% 选择一个示例试次进行PCA
trialForPCA = 1;
m2_activity_pca = squeeze(h_m2_all_test(trialForPCA, :, :));

% 执行PCA
[coeff, score, ~] = pca(m2_activity_pca);

% 绘制轨迹
scatter3(score(:,1), score(:,2), score(:,3), 30, time, 'filled');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
title(sprintf('M2状态空间轨迹 (试次%d)', trialForPCA));
colorbar;
colormap('jet');
grid on;
view(30, 20);

% 6. 模拟"M2失活"实验
subplot(3, 4, [11, 12]);

% 复制参数并减弱M2->下游连接
params_inactivated = params;
inactivation_factor = 0.1;
params_inactivated.W_m2_to_down = params.W_m2_to_down * inactivation_factor;

% 使用失活后的参数运行前向传播
[~, ~, ~, ~, outputs_inactivated] = forwardBackward(params_inactivated, testInputs, testTargets, baselineLength);

% 比较输出
plot(time, squeeze(testOutputs(trialIdx, :, 1)), 'b-', 'LineWidth', 2, 'DisplayName', '正常');
hold on;
plot(time, squeeze(outputs_inactivated(trialIdx, :, 1)), 'r--', 'LineWidth', 2, 'DisplayName', 'M2失活');
plot(time, squeeze(testTargets(trialIdx, :, 1)), 'k:', 'LineWidth', 1.5, 'DisplayName', '目标');

% 标记基线期
ylim = get(gca, 'YLim');
fill([1, baselineLength, baselineLength, 1], [ylim(1), ylim(1), ylim(2), ylim(2)], ...
    [0.9, 0.9, 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlabel('时间步'); ylabel('运动输出');
title('模拟M2失活实验 (无外部时间信号)');
legend('Location', 'best');
grid on;
hold off;

%% 8. 高级分析：检查内部计时机制
fprintf('\n分析内部计时机制...\n');

figure('Position', [100, 100, 1200, 500]);

% 计算基线期平均活动
baseline_mean = mean(h_m2_all_test(:, 1:baselineLength, :), 2);

% 计算相对于基线的变化
h_m2_delta = zeros(size(h_m2_all_test));
for b = 1:testBatchSize
    h_m2_delta(b, :, :) = squeeze(h_m2_all_test(b, :, :)) - squeeze(baseline_mean(b, 1, :))';
end

% 分析两种任务的内部动态
if ~isempty(cuedTrials) && ~isempty(selfTrials)
    % 线索触发任务：平均相对活动
    m2_delta_cued = mean(mean(h_m2_delta(cuedTrials, :, :), 3), 1);

    % 自主计时任务：平均相对活动
    m2_delta_self = mean(mean(h_m2_delta(selfTrials, :, :), 3), 1);

    % 绘制相对活动
    subplot(1, 3, 1);
    plot(time, m2_delta_cued, 'b-', 'LineWidth', 2, 'DisplayName', '线索触发');
    hold on;
    plot(time, m2_delta_self, 'r-', 'LineWidth', 2, 'DisplayName', '自主计时');

    % 标记基线期
    ylim = get(gca, 'YLim');
    fill([1, baselineLength, baselineLength, 1], [ylim(1), ylim(1), ylim(2), ylim(2)], ...
        [0.9, 0.9, 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    xlabel('时间步'); ylabel('相对于基线的活动变化');
    title('M2内部动态 (无外部时间信号)');
    legend('Location', 'best');
    grid on;
    hold off;

    % 分析动态模式
    subplot(1, 3, 2);

    % 计算任务期的关键特征
    task_time = (baselineLength+1):totalLength;

    % 线索触发：平台稳定性
    cued_end = m2_delta_cued(totalLength-19:totalLength);
    cued_stability = std(cued_end);

    % 自主计时：递增特征
    self_task = m2_delta_self(task_time);
    p = polyfit(task_time, self_task, 1);
    self_slope = p(1);

    % 计算达到平台的时间（线索触发）
    cued_task = m2_delta_cued(task_time);
    plateau_threshold = 0.8 * max(cued_task);
    plateau_time = find(cued_task >= plateau_threshold, 1);
    if isempty(plateau_time)
        plateau_time = length(task_time);
    end

    bar([1, 2, 3], [cued_stability, abs(self_slope), plateau_time], 'FaceColor', [0.7, 0.7, 0.9]);
    set(gca, 'XTickLabel', {'平台稳定性', '递增斜率', '达到平台时间'});
    ylabel('指标值');
    title('内部计时特征');
    grid on;

    % 3. 时间解码分析
    subplot(1, 3, 3);

    % 使用线性回归解码时间
    X = squeeze(mean(h_m2_delta(cuedTrials, baselineLength+1:end, :), 1)); % [taskLength, n_m2]
    Y = (1:taskLength)'; % 真实时间

    % 训练线性解码器
    beta = (X' * X) \ (X' * Y);
    Y_pred = X * beta;

    plot(1:taskLength, Y, 'k-', 'LineWidth', 2, 'DisplayName', '真实时间');
    hold on;
    plot(1:taskLength, Y_pred, 'b--', 'LineWidth', 1.5, 'DisplayName', '解码时间');
    xlabel('任务期时间步'); ylabel('时间估计');
    title('从M2活动解码时间');
    legend('Location', 'best');
    grid on;

    % 计算解码准确度
    r2 = 1 - sum((Y - Y_pred).^2) / sum((Y - mean(Y)).^2);
    fprintf('时间解码准确度 (R²): %.4f\n', r2);
    fprintf('内部计时特征:\n');
    fprintf('  线索触发 - 平台稳定性(标准差): %.4f\n', cued_stability);
    fprintf('  自主计时 - 递增斜率: %.4f\n', self_slope);
    fprintf('  线索触发 - 达到平台时间: %d 步\n', plateau_time);
end

% %% 9. 保存结果
% save('m2_rnn_no_timesignal.mat', 'params', 'lossHistory', 'testLoss', 'baselineLength', 'totalLength', 'inputSize');
% fprintf('\n模型已保存到 m2_rnn_no_timesignal.mat\n');
% fprintf('无外部时间信号模型运行完成！\n');

%% 2. 数据生成函数 (无外部时间信号)
function [inputs, targets, trialInfo] = generateBatch(batchSize, totalLength, baselineLength, cueProb)
% 生成一批训练数据（无外部时间信号）
% 输出:
%   inputs: [batchSize, totalLength, 2] (仅任务指令)
%   targets: [batchSize, totalLength, 1]
%   trialInfo: 结构数组，包含每个试次的信息

inputSize = 2;  % 仅任务指令
outputSize = 1;
taskLength = totalLength - baselineLength;

inputs = zeros(batchSize, totalLength, inputSize);
targets = zeros(batchSize, totalLength, outputSize);
trialInfo = struct('isCued', cell(batchSize, 1), 'forePeriod', [], 'pulseTime', []);

for b = 1:batchSize
    % 随机决定任务类型
    isCued = rand() < cueProb;
    trialInfo(b).isCued = isCued;

    % 随机生成任务期内的准备期长度
%     forePeriodInTask = randi([20, taskLength-20]);
    forePeriods = [40 80 120];
    forePeriodInTask = forePeriods(randi(3));
    trialInfo(b).forePeriod = forePeriodInTask;

    % 计算在整个试次中的绝对运动触发时刻
    pulseTimeAbs = baselineLength + forePeriodInTask;
    trialInfo(b).pulseTime = pulseTimeAbs;

    % ==== 设置基线期 (t = 1:baselineLength) ====
    for t = 1:baselineLength
        % 基线期: 无任务指令
        inputs(b, t, :) = [0, 0];    % 仅任务指令，无时间信号
        targets(b, t, :) = 0;        % 基线期无运动输出
    end

    % ==== 设置任务期 ====
    transitionSteps = 15; % 过渡期长度，建议10-20步

    for t_rel = 1:taskLength
        t_abs = baselineLength + t_rel;

        % 计算平滑过渡系数
        if t_rel <= transitionSteps
            % 在过渡期内使用平滑函数
            ramp_factor = sigmoid_ramp(t_rel, transitionSteps); % 使用方案2
        else
            % 过渡期结束后，指令保持全强度
            ramp_factor = 1;
        end

        if isCued
            % 线索触发任务: [ramp_factor, 0]
            inputs(b, t_abs, :) = [ramp_factor, 0];
        else
            % 自主计时任务: [0, ramp_factor]
            inputs(b, t_abs, :) = [0, ramp_factor];
        end
    end

    % ==== 设置运动输出目标 (高斯脉冲) ====
    pulseWidth = 3;
    for t_abs = 1:totalLength
        targets(b, t_abs, 1) = exp(-((t_abs - pulseTimeAbs)^2) / (2 * pulseWidth^2));
    end
end

inputs = single(inputs);
targets = single(targets);
end

%% 平滑过渡函数定义
% 方案1: 线性渐变 (最简)
function s = linear_ramp(t, totalSteps)
    % t: 当前步 (1到totalSteps)
    % totalSteps: 总过渡步数
    s = min(max(t / totalSteps, 0), 1);
end

% 方案2: S型（Sigmoid）渐变 (最符合神经激活特性)
function s = sigmoid_ramp(t, totalSteps)
    % 使用sigmoid函数实现平滑启停
    slope = 8; % 控制过渡陡度，越大越接近阶跃
    s = 1 ./ (1 + exp(-slope * (2*t/totalSteps - 1)));
end

% 方案3: 双指数渐变 (模拟快速上升与缓慢调整)
function s = double_exp_ramp(t, totalSteps)
    % 快速上升与缓慢稳定的组合
    rise_factor = 0.7;
    adjust_factor = 0.3;
    s = (1 - exp(-t/(totalSteps*rise_factor))) * (1 + adjust_factor*exp(-t/(totalSteps*0.5)));
    s = min(s, 1); % 限制最大值
end

%% 4. 前向传播与反向传播函数
function [loss, grads, h_m2_all, h_down_all, outputs] = forwardBackward(params, inputs, targets, baselineLength)
% 解包参数
[W_m2_rec, W_m2_in, b_m2, W_m2_to_down, W_down_rec, b_down, W_out, b_out] = ...
    deal(params.W_m2_rec, params.W_m2_in, params.b_m2, params.W_m2_to_down, ...
    params.W_down_rec, params.b_down, params.W_out, params.b_out);

[batchSize, totalLength, inputDim] = size(inputs);
n_m2 = size(W_m2_rec, 1);
n_down = size(W_down_rec, 1);
inputSize = size(W_m2_in, 2);
outputSize = size(W_out, 1);

% 验证输入维度
if inputDim ~= inputSize
    error('输入维度不匹配: inputs维度=%d, W_m2_in期望=%d', inputDim, inputSize);
end

% 初始化隐藏状态
h_m2 = zeros(n_m2, batchSize);
h_down = zeros(n_down, batchSize);

% 存储所有时间步的状态
h_m2_all = zeros(batchSize, totalLength, n_m2);
h_down_all = zeros(batchSize, totalLength, n_down);
outputs = zeros(batchSize, totalLength, outputSize);

% ===== 前向传播 =====
loss = 0;
taskStepsCount = 0;

for t = 1:totalLength
    % 获取当前时间步输入 [inputSize, batchSize]
    input_t = reshape(inputs(:, t, :), [], batchSize);
    input_t = reshape(input_t, [inputSize, batchSize]);

    % M2模块更新
    h_m2 = tanh(W_m2_rec * h_m2 + W_m2_in * input_t + b_m2);

    % 下游模块更新
    h_down = tanh(W_down_rec * h_down + W_m2_to_down * h_m2 + b_down);

    % 生成输出
    y_t = W_out * h_down + b_out;
    outputs(:, t, :) = reshape(y_t, [batchSize, 1, outputSize]);

    % 存储状态
    h_m2_all(:, t, :) = reshape(h_m2, [batchSize, 1, n_m2]);
    h_down_all(:, t, :) = reshape(h_down, [batchSize, 1, n_down]);

    % 只在任务期计算损失
    if t > baselineLength
        target_t = reshape(targets(:, t, :), [], batchSize);
        target_t = reshape(target_t, [outputSize, batchSize]);
        loss = loss + sum((y_t - target_t).^2, 'all');
        taskStepsCount = taskStepsCount + 1;
    end
end

% 平均损失 (只基于任务期)
if taskStepsCount > 0
    loss = loss / (batchSize * taskStepsCount);
else
    loss = 0;
end

% ===== 反向传播 (BPTT) =====
grads = struct('W_m2_rec', zeros(size(W_m2_rec)), 'W_m2_in', zeros(size(W_m2_in)), 'b_m2', zeros(size(b_m2)), ...
    'W_m2_to_down', zeros(size(W_m2_to_down)), ...
    'W_down_rec', zeros(size(W_down_rec)), 'b_down', zeros(size(b_down)), ...
    'W_out', zeros(size(W_out)), 'b_out', zeros(size(b_out)));

dh_m2_next = zeros(size(h_m2));
dh_down_next = zeros(size(h_down));

% 只在任务期回传梯度
for t = totalLength:-1:(baselineLength+1)
    % 获取当前时间步的状态
    h_m2_t = reshape(h_m2_all(:, t, :), [], batchSize);
    h_m2_t = reshape(h_m2_t, [n_m2, batchSize]);

    h_down_t = reshape(h_down_all(:, t, :), [], batchSize);
    h_down_t = reshape(h_down_t, [n_down, batchSize]);

    % 获取前一时间步状态
    if t > 1
        h_m2_prev = reshape(h_m2_all(:, t-1, :), [], batchSize);
        h_m2_prev = reshape(h_m2_prev, [n_m2, batchSize]);

        h_down_prev = reshape(h_down_all(:, t-1, :), [], batchSize);
        h_down_prev = reshape(h_down_prev, [n_down, batchSize]);
    else
        h_m2_prev = zeros(n_m2, batchSize);
        h_down_prev = zeros(n_down, batchSize);
    end

    % 获取输入
    input_t = reshape(inputs(:, t, :), [], batchSize);
    input_t = reshape(input_t, [inputSize, batchSize]);

    % 获取输出和目标
    y_t = reshape(outputs(:, t, :), [], batchSize);
    y_t = reshape(y_t, [outputSize, batchSize]);

    target_t = reshape(targets(:, t, :), [], batchSize);
    target_t = reshape(target_t, [outputSize, batchSize]);

    % 输出层梯度
    dL_dy = 2 * (y_t - target_t) / (batchSize * taskStepsCount);

    grads.W_out = grads.W_out + dL_dy * h_down_t';
    grads.b_out = grads.b_out + sum(dL_dy, 2);

    % 下游模块梯度
    dL_dh_down = W_out' * dL_dy + dh_down_next;
    dtanh_down = 1 - h_down_t.^2;
    delta_down = dL_dh_down .* dtanh_down;

    grads.W_down_rec = grads.W_down_rec + delta_down * h_down_prev';
    grads.b_down = grads.b_down + sum(delta_down, 2);
    grads.W_m2_to_down = grads.W_m2_to_down + delta_down * h_m2_t';

    % M2模块梯度
    dL_dh_m2 = W_m2_to_down' * delta_down + dh_m2_next;
    dtanh_m2 = 1 - h_m2_t.^2;
    delta_m2 = dL_dh_m2 .* dtanh_m2;

    grads.W_m2_rec = grads.W_m2_rec + delta_m2 * h_m2_prev';
    grads.W_m2_in = grads.W_m2_in + delta_m2 * input_t';
    grads.b_m2 = grads.b_m2 + sum(delta_m2, 2);

    % 为上一时间步准备梯度
    dh_down_next = W_down_rec' * delta_down;
    dh_m2_next = W_m2_rec' * delta_m2;
end
end
