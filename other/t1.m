% %--------parameter1--------
% l_max = 600;
% l_min = 300;
% vA1 = [100 0 0];
% vA2 = [100*cos(2*pi/3) 100*sin(2*pi/3) 0];
% vA3 = [100*cos(4*pi/3) 100*sin(4*pi/3) 0];
% va1 = [50 0 0];
% va2 = [50*cos(2*pi/3) 50*sin(2*pi/3) 0];
% va3 = [50*cos(4*pi/3) 50*sin(4*pi/3) 0];
% VA = [vA1; vA2; vA3];
% Va = [va1; va2; va3];
% % ------search space-------
% seq_x = -600 : 10 : 600;
% seq_y = -600 : 10 : 600;
% seq_z = 50 : 5 : 650;
% % -----end-parameter1------

% -------parameter2------
l_max = 504.5;
l_min = 454.5;
xa=[92.58 132.58 40 -40 -132.58 -92.58];
ya=[99.64 30.36 -130 -130 30.36 99.64];
za=23.1*ones(1,6);
xb=[30 78.22 48.22 -48.22 -78.22 -30];
yb=[73 -10.52 -62.48 -62.48 -10.52 73];
zb=-37.1*ones(1,6);
VA = [xa', ya', za'];
Va = [xb', yb', zb'];
% -------search space------
seq_x = -300 : 5 : 300;
seq_y = -300 : 5 : 300;
seq_z = 400 : 2 : 650;
% -----end-parameter2------

% %--------parameter3--------
% l_max = 1000;
% l_min = 650;
% R1 = 800;
% R2 = 600;
% H = 20;
% r1 = 100;
% r2 = 80;
% h = 100;
% B1 = [R1*cos(pi/2);   R1*sin(pi/2);   0];
% B2 = [R1*cos(7*pi/6); R1*sin(7*pi/6); 0];
% B3 = [R1*cos(-pi/6);  R1*sin(-pi/6);  0];
% B4 = [R2*cos(pi/6);   R2*sin(pi/6);   H];
% B5 = [R2*cos(5*pi/6); R2*sin(5*pi/6); H];
% P1_m = [r1*cos(pi/2);   r1*sin(pi/2);   0];
% P2_m = [r1*cos(7*pi/6); r1*sin(7*pi/6); 0];
% P3_m = [r1*cos(-pi/6);  r1*sin(-pi/6);  0];
% P4_m = [r2*cos(pi/6);   r2*sin(pi/6);   h];
% P5_m = [r2*cos(5*pi/6); r2*sin(5*pi/6); h];
% VA = [B1 B2 B3 B4 B5]';
% Va = [P1_m P2_m P3_m P4_m P5_m]';
% % ------search space-------
% seq_x = -600 : 10 : 600;
% seq_y = -600 : 10 : 600;
% seq_z = -800 : 10 : -200;
% % -----end-parameter3------


% assistant parameter
work_space = [0;0;0];


% length transform
for ix = 1 : length(seq_x)
    for iy = 1 : length(seq_y)
        z_min_point = zeros(3,1);
        z_max_point = zeros(3,1);
        for iz = 1 : length(seq_z)
            vt = [seq_x(ix) seq_y(iy) seq_z(iz)];  % 搜索的目标点
            pos_flag = 0;  % 位置可达标志位
            
            for j = 1 : length(Va(:,1))
                vAa = vt + Va(j, :) - VA(j, :);
                len_vAa = norm(vAa);

                if (len_vAa >= l_min)&&(len_vAa <= l_max)
                    pos_flag = pos_flag + 1;
                end

                if pos_flag == length(Va(:,1))
                    if z_min_point == zeros(3,1)
                        z_min_point = vt';
                        z_max_point = vt';
                    else
                        z_max_point = vt';
                   end
                end

                
            end        
        end
        if (z_min_point(3) ~= 0) && (z_max_point(3) ~= 0)
            work_space = [work_space z_min_point z_max_point];
        end
        
    end
end

figure('Color', [1 1 1])
plot3(work_space(1,2:end), work_space(2,2:end), work_space(3,2:end), '.', 'Color', '#87CEFA');
% plot3(vt(1), vt(2), vt(3), '.', 'Color', '#87CEFA');
grid on
fprintf('done\n');









