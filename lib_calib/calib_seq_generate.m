function [Pos_ref_seq, Pos_valid_seq] = calib_seq_generate(unit_para)
%CALIB_SEQ_GENERATE 生成标定位姿参考序列和验证序列
%   输入：
%     unit_para — 单位转换因子（m/单位），如 basic_read 返回的 unit_para
%   输出：
%     Pos_ref_seq   — 5×nref 参考位姿序列 [x; y; z; phi; theta]（rad）
%     Pos_valid_seq — 5×nvalid 验证位姿序列 [x; y; z; phi; theta]（rad）
%
%   参考序列通过网格采样生成，验证序列通过随机采样生成。
%   跳过了 theta=0 且 phi~=0 的冗余组合（姿态奇异性区域）。

%% 参考序列：网格采样
x_seq = [-100 0 100] * unit_para;
y_seq = [-100 0 100] * unit_para;
z_seq = [-900 -1000] * unit_para;
phi_seq = deg2rad(-120 : 120 : 120);
theta_seq = deg2rad([0 5 10]);

Pos_ref_seq = zeros(5, 3*3*2*3*3);
idx = 0;
for ix = 1:length(x_seq)
    for iy = 1:length(y_seq)
        for iz = 1:length(z_seq)
            for iphi = 1:length(phi_seq)
                for itheta = 1:length(theta_seq)
                    if theta_seq(itheta) == 0 && phi_seq(iphi) ~= 0
                        continue;
                    end
                    idx = idx + 1;
                    Pos_ref_seq(:, idx) = [x_seq(ix); y_seq(iy); z_seq(iz); phi_seq(iphi); theta_seq(itheta)];
                end
            end
        end
    end
end
Pos_ref_seq = Pos_ref_seq(:, 1:idx);

%% 验证序列：随机采样
x_v = 200 * unit_para * (rand(3, 1) - 0.5);
y_v = 200 * unit_para * (rand(3, 1) - 0.5);
z_v = (-900 + 100 * rand(3, 1)) * unit_para;
phi_v = deg2rad(300 * (rand(3, 1) - 0.5));
theta_v = deg2rad(10 * rand(3, 1));

Pos_valid_seq = zeros(5, 3^5);
idx = 0;
for ix = 1:3, for iy = 1:3, for iz = 1:3, for iphi = 1:3, for itheta = 1:3
    if theta_v(itheta) == 0 && phi_v(iphi) ~= 0
        continue;
    end
    idx = idx + 1;
    Pos_valid_seq(:, idx) = [x_v(ix); y_v(iy); z_v(iz); phi_v(iphi); theta_v(itheta)];
end, end, end, end, end
Pos_valid_seq = Pos_valid_seq(:, 1:idx);

end
