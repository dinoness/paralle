clear
addpath(genpath('./lib'));
%% 参数集
%--------parameter3--------
unit_para = 1;  % 0.001表示m，1表示mm

T = readtable('parameters.xlsx', 'Range', 'A2:B12');
paras = table2array(T(:, 2));
l_max = paras(1)*unit_para;
l_min = paras(2)*unit_para;  % 670
R1 = paras(3)*unit_para;  % 550
R2 = paras(4)*unit_para;  % 500
H = paras(5)*unit_para;  % 0
r1 = paras(6)*unit_para;  % 100
r2 = paras(7)*unit_para;  % 80
h = paras(8)*unit_para;  % 10

limb_dir = [pi/2; 7*pi/6; -pi/6; pi/6; 5*pi/6];
B1 = [R1*cos(pi/2);   R1*sin(pi/2);   0];
B2 = [R1*cos(7*pi/6); R1*sin(7*pi/6); 0];
B3 = [R1*cos(-pi/6);  R1*sin(-pi/6);  0];
B4 = [R2*cos(pi/6);   R2*sin(pi/6);   H];
B5 = [R2*cos(5*pi/6); R2*sin(5*pi/6); H];
B = [B1 B2 B3 B4 B5];

% move plant parameter
P1_m = [r1*cos(pi/2);   r1*sin(pi/2);   0];
P2_m = [r1*cos(7*pi/6); r1*sin(7*pi/6); 0];
P3_m = [r1*cos(-pi/6);  r1*sin(-pi/6);  0];
P4_m = [r2*cos(pi/6);   r2*sin(pi/6);   h];
P5_m = [r2*cos(5*pi/6); r2*sin(5*pi/6); h];
P_m = [P1_m P2_m P3_m P4_m P5_m];
P_v = zeros(3, 5);  % 只变换了方向，没变换起点
P = zeros(3, 5);    % 末端点坐标

% -----end-parameter3------

Pos_ref_seq = [0*unit_para;0*unit_para;-600*unit_para;0;0];  % line=5 colum=n  角度的单位是° **一列为一组**
T_ref = pos2trans(Pos_ref_seq(:, 1), B);
l0 = 600;
l0_seq = [l0;l0;l0;l0;l0];
joint_u_angle_tilt = 155 / 180 * pi;
p_seq = parameterize(limb_dir, B, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

joint_q0 = keni_sol_inverse(T_ref, B, l0_seq, P_m, p_seq);


B_delta = B;
B_delta(3,3) = B_delta(3,3) + 3;
B_delta(2,3) = B_delta(2,3) - 3;
p_seq2 = parameterize(limb_dir, B_delta, r1, r2, l0_seq, P_m, joint_u_angle_tilt);

T_delta = keni_sol_forward(joint_q0, p_seq2, 1e-8)


% 作图
draw_plate = [32 62 111 111 79 49 -49 -79 -111 -111 -62 -32 32;
            110 93 7 -27 -83 -100 -100 -83 -27 7 93 110 110;
            -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5 -31.5];
darw_bar = [0 0; 0 0; -31.5 -100];

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
% plot3(P(1, :), P(2, :), -1*P(3, :), 'o', 'Color', '#4682B4');
B_plot = [B(:,1) B(:,5) B(:,2) B(:,3) B(:,4) B(:,1)];
% P_plot = [P(:,1) P(:,5) P(:,2) P(:,3) P(:,4) P(:,1)];
plot3(B_plot(1, :), B_plot(2, :), -1*B_plot(3, :), '-', 'Color', '#4682B4');
% plot3(P_plot(1, :), P_plot(2, :), -1*P_plot(3, :), '-', 'Color', '#4682B4');
plot3(draw_p_origin(1, :), draw_p_origin(2, :), -1*draw_p_origin(3, :), '-', 'Color', '#4682B4');
plot3(draw_b_origin(1, :), draw_b_origin(2, :), -1*draw_b_origin(3, :), '-', 'Color', '#4682B4');

% for i = 1 : 5
%     plot3([B(1, i) P(1, i)], [B(2, i) P(2, i)], -1*[B(3, i) P(3, i)], '-', 'Color', '#4682B4');
% end
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
% plot3(P_delta(1, :), P_delta(2, :), -1*P_delta(3, :), 'o', 'Color', '#FF7F50');
B_plot2 = [B_delta(:,1) B_delta(:,5) B_delta(:,2) B_delta(:,3) B_delta(:,4) B_delta(:,1)];
% P_plot2 = [P_delta(:,1) P_delta(:,5) P_delta(:,2) P_delta(:,3) P_delta(:,4) P_delta(:,1)];
plot3(B_plot2(1, :), B_plot2(2, :), -1*B_plot2(3, :), '-', 'Color', '#FF7F50');
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






