% 指定结构误差参数，查看位姿变化
clear
flag_plot = 1;

path_add();
basic_paras = basic_read('parameters.xlsx');
l_max = basic_paras.l_max;
l_min = basic_paras.l_min;
R1 = basic_paras.R1;
R2 = basic_paras.R2;
H = basic_paras.H;
r1 = basic_paras.r1;
r2 = basic_paras.r2;
h = basic_paras.h;
L_tool = basic_paras.L_tool;
limb_dir = basic_paras.limb_dir;
B = basic_paras.B;
P_m = basic_paras.P_m;
l0_seq = basic_paras.l0_seq;


unit_para = 0.001;  % 0.001表示m，1表示mm
% -----end-parameter3------



P_v = zeros(3, 5);  % 只变换了方向，没变换起点
P = zeros(3, 5);    % 末端点坐标
Pos_ref_seq = [0*unit_para;0*unit_para;-800*unit_para;0;0];  % line=5 colum=n  角度的单位是° **一列为一组**
T_ref = pos2trans(Pos_ref_seq(:, 1), B);
joint_u_angle_tilt = 155 / 180 * pi;
p_seq = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

joint_q0 = keni_sol_inverse(T_ref, B, l0_seq, P_m, p_seq);



% B_delta(3,3) = B_delta(3,3) - 0.01;
A = [-1 0 0;0 0 1;0 1 0];
% 施加5Nm扭矩，钢架形变
% B_delta = A \ [0.15531 0.96248  0.95566  0.41285 0.41441;
%                -2.8757e-002 0.66129 -0.68776 -0.52018 0.51169;
%                -2.6366e-002 0.47564  -0.46752  -0.42598 0.42688]*1e-3;

% 施加5Nm扭矩，钢架+型钢强化形变
% B_delta = A \ [1.8657e-002 0.12832  0.13184  5.3306e-002 5.2661e-002;
%                -1.1315e-004 0.17818 -0.20102 -6.8149e-002 6.145e-002;
%                2.6872e-003 6.0572e-002  -5.822e-002  -4.9764e-002 5.2126e-002]*1e-3;

% v3大理石架施加5Nm扭矩+30N力，形变量5.06μm
% B_delta = A \ [-5.8614e-005 6.9772e-004  3.9442e-005 -1.891e-004  -3.3246e-005;
%                -5.949e-004 -1.041e-003 -8.5822e-004 -1.0678e-003 -9.0932e-004;
%                 -1.3539e-004 -2.1966e-004  -8.6195e-004 -3.7212e-004  3.0401e-004]*1e-3;


% B_delta = A \ [-1.5098e-005 -4.3266e-004  4.9412e-004  1.0722e-005 -5.5064e-005;
%                -5.7798e-005 -3.1701e-005 -5.2699e-005 -7.1807e-006 -2.9126e-005;
%                 4.2404e-005 -4.1589e-005  4.8952e-005  9.2001e-005 -1.0792e-005]*1e-3;
% B_delta = A \ [-1.4214e-005 -4.6059e-004  5.4223e-004  2.1144e-005 -6.6548e-005;
%                -5.2966e-005 -3.5028e-005 -5.8474e-005 -5.3564e-006 -3.0587e-005;
%                 4.904e-005 -4.4051e-005  5.8018e-005  1.0583e-004 -9.1656e-006]*1e-3;
B_delta = zeros(3, 5);
B_delta = B + B_delta;
% B_delta(3,2) = B_delta(3,2) + 0.015/1000;
% B_delta(2,2) = B_delta(2,2) + 0.015/1000;
% B_delta(3,1) = B_delta(3,1) - 0.015/1000;

p_seq2 = parameterize(limb_dir, B_delta, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

T_delta = keni_sol_forward(joint_q0, p_seq2, 1e-8);
format long
err = T_ref\T_delta;
disp(norm(err(1:3, 4)))


if flag_plot == 1
    % 作图
    draw_plate = [32 62 111 111 79 49 -49 -79 -111 -111 -62 -32 32;
                110 93 7 -27 -83 -100 -100 -83 -27 7 93 110 110;
                -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5] * unit_para;
    darw_bar = [0 0; 0 0; -31.5 -100] * unit_para;

    draw_p_origin = zeros(3,13);
    draw_b_origin = zeros(3, 2);

    fig = figure('Color', [1 1 1]);
    % 机构简图
    for i = 1 : 5
        P(:, i) = T_ref(1:3, 1:3) * P_m(:, i) + T_ref(1:3, 4);
    end
    for i = 1 : 13
        draw_p_origin(:, i) = T_ref(1:3, 1:3) * draw_plate(:, i) + T_ref(1:3, 4);
    end
    for i = 1 : 2
        draw_b_origin(:, i) =  T_ref(1:3, 1:3) * darw_bar(:, i) + T_ref(1:3, 4);
    end
    plot3(B(1, :), B(2, :), -1*B(3, :), 'o', 'Color', '#4682B4');
    hold on
    plot3(P(1, :), P(2, :), -1*P(3, :), 'o', 'Color', '#4682B4');
    B_plot = [B(:,1) B(:,5) B(:,2) B(:,3) B(:,4) B(:,1)];
    P_plot = [P(:,1) P(:,5) P(:,2) P(:,3) P(:,4) P(:,1)];

    plot3(B_plot(1, :), B_plot(2, :), -1*B_plot(3, :), '-', 'Color', '#4682B4');
    plot3(P_plot(1, :), P_plot(2, :), -1*P_plot(3, :), '-', 'Color', '#4682B4');
    plot3(draw_p_origin(1, :), draw_p_origin(2, :), -1*draw_p_origin(3, :), '-', 'Color', '#4682B4','LineWidth',1.5);
    hold on
    plot3(draw_b_origin(1, :), draw_b_origin(2, :), -1*draw_b_origin(3, :), '-', 'Color', '#4682B4','LineWidth',1.5);

    for i = 1 : 5
        plot3([B(1, i) P(1, i)], [B(2, i) P(2, i)], -1*[B(3, i) P(3, i)], '-', 'Color', '#4682B4');
    end
    % ============
    P_delta = zeros(3,5);
    for i = 1 : 5
        P_delta(:, i) = T_delta(1:3, 1:3) * P_m(:, i) + T_delta(1:3, 4);
    end

    for i = 1 : 13
        draw_p_origin(:, i) = T_delta(1:3, 1:3) * draw_plate(:, i) + T_delta(1:3, 4);
    end
    for i = 1 : 2
        draw_b_origin(:, i) =  T_delta(1:3, 1:3) * darw_bar(:, i) + T_delta(1:3, 4);
    end

    plot3(B_delta(1, :), B_delta(2, :), -1*B_delta(3, :), 'o', 'Color', '#FF7F50')
    plot3(P_delta(1, :), P_delta(2, :), -1*P_delta(3, :), 'o', 'Color', '#FF7F50');
    B_plot2 = [B_delta(:,1) B_delta(:,5) B_delta(:,2) B_delta(:,3) B_delta(:,4) B_delta(:,1)];
    P_plot2 = [P_delta(:,1) P_delta(:,5) P_delta(:,2) P_delta(:,3) P_delta(:,4) P_delta(:,1)];
    plot3(B_plot2(1, :), B_plot2(2, :), -1*B_plot2(3, :), '-', 'Color', '#FF7F50');
    plot3(P_plot2(1, :), P_plot2(2, :), -1*P_plot2(3, :), '-', 'Color', '#FF7F50');
    for i = 1 : 5
        plot3([B_delta(1, i) P_delta(1, i)], [B_delta(2, i) P_delta(2, i)], -1*[B_delta(3, i) P_delta(3, i)], '-', 'Color', '#FF7F50');
    end

    plot3(draw_p_origin(1, :), draw_p_origin(2, :), -1*draw_p_origin(3, :), '-', 'Color', '#FF7F50','LineWidth',1.5);
    plot3(draw_b_origin(1, :), draw_b_origin(2, :), -1*draw_b_origin(3, :), '-', 'Color', '#FF7F50','LineWidth',1.5);

    % fig = figure('Color', [1 1 1]);
    % % 机构简图
    % for i = 1 : 5
    %     P(:, i) = T_ref(1:3, 1:3) * P_m(:, i) + T_ref(1:3, 4);
    % end
    % for i = 1 : 13
    %     draw_p_origin(:, i) = T_ref(1:3, 1:3) * draw_plate(:, i) + T_ref(1:3, 4);
    % end
    % for i = 1 : 2
    %     draw_b_origin(:, i) =  T_ref(1:3, 1:3) * darw_bar(:, i) + T_ref(1:3, 4);
    % end
    % plot3(B(1, :), B(2, :), -1*B(3, :), 'o', 'Color', '#4682B4');
    % hold on
    % plot3(P(1, :), P(2, :), -1*P(3, :), 'o', 'Color', '#4682B4');
    % B_plot = [B(:,1) B(:,5) B(:,2) B(:,3) B(:,4) B(:,1)];
    % P_plot = [P(:,1) P(:,5) P(:,2) P(:,3) P(:,4) P(:,1)];

    % plot3(B_plot(1, :), B_plot(2, :), -1*B_plot(3, :), '-', 'Color', '#4682B4');
    % plot3(P_plot(1, :), P_plot(2, :), -1*P_plot(3, :), '-', 'Color', '#4682B4');
    % plot3(draw_p_origin(1, :), draw_p_origin(2, :), -1*draw_p_origin(3, :), '-', 'Color', '#4682B4');
    % hold on
    % plot3(draw_b_origin(1, :), draw_b_origin(2, :), -1*draw_b_origin(3, :), '-', 'Color', '#4682B4');

    % for i = 1 : 5
    %     plot3([B(1, i) P(1, i)], [B(2, i) P(2, i)], -1*[B(3, i) P(3, i)], '-', 'Color', '#4682B4');
    % end
    % % ============
    % P_delta = zeros(3,5);
    % for i = 1 : 5
    %     P_delta(:, i) = T_delta(1:3, 1:3) * P_m(:, i) + T_delta(1:3, 4);
    % end

    % for i = 1 : 13
    %     draw_p_origin(:, i) = T_delta(1:3, 1:3) * draw_plate(:, i) + T_delta(1:3, 4);
    % end
    % for i = 1 : 2
    %     draw_b_origin(:, i) =  T_delta(1:3, 1:3) * darw_bar(:, i) + T_delta(1:3, 4);
    % end

    % plot3(B_delta(1, :), B_delta(2, :), -1*B_delta(3, :), 'o', 'Color', '#FF7F50')
    % plot3(P_delta(1, :), P_delta(2, :), -1*P_delta(3, :), 'o', 'Color', '#FF7F50');
    % B_plot2 = [B_delta(:,1) B_delta(:,5) B_delta(:,2) B_delta(:,3) B_delta(:,4) B_delta(:,1)];
    % P_plot2 = [P_delta(:,1) P_delta(:,5) P_delta(:,2) P_delta(:,3) P_delta(:,4) P_delta(:,1)];
    % plot3(B_plot2(1, :), B_plot2(2, :), -1*B_plot2(3, :), '-', 'Color', '#FF7F50');
    % plot3(P_plot2(1, :), P_plot2(2, :), -1*P_plot2(3, :), '-', 'Color', '#FF7F50');
    % for i = 1 : 5
    %     plot3([B_delta(1, i) P_delta(1, i)], [B_delta(2, i) P_delta(2, i)], -1*[B_delta(3, i) P_delta(3, i)], '-', 'Color', '#FF7F50');
    % end

    plot3(draw_p_origin(1, :), draw_p_origin(2, :), -1*draw_p_origin(3, :), '-', 'Color', '#FF7F50');
    plot3(draw_b_origin(1, :), draw_b_origin(2, :), -1*draw_b_origin(3, :), '-', 'Color', '#FF7F50');
    % #4682B4 #32CD32 #FF7F50
    grid on
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')


end



