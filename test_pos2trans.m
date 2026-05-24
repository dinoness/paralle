path_add()
basic_paras = basic_read('parameters.xlsx', 'column', 'B', 'unit', 'mm');
B = basic_paras.B;

% Test 1: default behavior (expects deg)
Pos_ref_deg = [0; 0; -800; 90; 90];  % phi=90, theta=90 in deg
T1 = pos2trans(Pos_ref_deg, B);
R1 = T1(1:3, 1:3);
z_axis1 = R1(:, 3);

% Test 2: with rad unit
Pos_ref_rad = [0; 0; -800; pi/2; pi/2];  % phi=90, theta=90 in rad
T2 = pos2trans(Pos_ref_rad, B, 'unit', 'rad');
R2 = T2(1:3, 1:3);
z_axis2 = R2(:, 3);

fprintf('Test 1 (deg input, default): z_axis = [%.4f, %.4f, %.4f]\n', z_axis1);
fprintf('Test 2 (rad input, unit=rad): z_axis = [%.4f, %.4f, %.4f]\n', z_axis2);
fprintf('Expected z_axis for theta=90, phi=90: [0, 1, 0] (approx)\n');

% Test 3: the bug - rad input but default unit
T3 = pos2trans(Pos_ref_rad, B);
R3 = T3(1:3, 1:3);
z_axis3 = R3(:, 3);
fprintf('Test 3 (BUG: rad input, default): z_axis = [%.4f, %.4f, %.4f]\n', z_axis3);
fprintf('This corresponds to theta=%.4f deg, not 90 deg!\n', norm(Pos_ref_rad(5))/pi*180);
