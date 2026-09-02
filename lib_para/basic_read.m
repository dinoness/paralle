function basic_paras = basic_read(file_name, varargin)
    range_start = 2;
    range_end = 30;
    column = 'B';
    unit_para = 0.001;  % 0.001表示m，1表示mm
    
    for k = 1:2:numel(varargin)
        switch lower(varargin{k})
            case 'column'
                column = varargin{k+1};
            case 'unit'
                switch (varargin{k+1})
                    case 'mm'
                        unit_para = 1;
                    case 'm'
                        unit_para = 0.001;
                    otherwise
                        error("Undefined unit: %s", varargin{k+1})
                end
            otherwise
                error('Unknown option: %s', varargin{k});
        end
    end

    % range = 'B2:B30';
    range = sprintf('%s%d:%s%d', column, range_start, column, range_end);

    T = readtable(file_name, 'Range', range);
    paras = table2array(T(:, 1));
    l_max = paras(1)*unit_para;
    l_min = paras(2)*unit_para;  % 670
    R1 = paras(3)*unit_para;  % 550
    R2 = paras(4)*unit_para;  % 500
    H = paras(5)*unit_para;  % 0
    r1 = paras(6)*unit_para;  % 100
    r2 = paras(7)*unit_para;  % 80
    h = paras(8)*unit_para;  % 10
    L_tool = paras(9)*unit_para;
    joint_u_angle_tilt = deg2rad(paras(10));
    % L_tool = 0;
    limb_dir = zeros(5, 2);
    l0_seq = zeros(5, 1);
    for i = 1 : 5
        limb_dir(i, 1) = deg2rad(paras(10+i));
        limb_dir(i, 2) = deg2rad(paras(10+5+i));
        l0_seq(i) = paras(10+5+5+i)*unit_para;
    end

    B1 = [R1*cos(limb_dir(1,1)); R1*sin(limb_dir(1,1)); 0];
    B2 = [R1*cos(limb_dir(2,1)); R1*sin(limb_dir(2,1)); 0];
    B3 = [R1*cos(limb_dir(3,1)); R1*sin(limb_dir(3,1)); 0];
    B4 = [R2*cos(limb_dir(4,1)); R2*sin(limb_dir(4,1)); H];
    B5 = [R2*cos(limb_dir(5,1)); R2*sin(limb_dir(5,1)); H];
    B = [B1 B2 B3 B4 B5];

    % move plant parameter
    P1_m = [r1*cos(limb_dir(1,2)); r1*sin(limb_dir(1,2)); L_tool];
    P2_m = [r1*cos(limb_dir(2,2)); r1*sin(limb_dir(2,2)); L_tool];
    P3_m = [r1*cos(limb_dir(3,2)); r1*sin(limb_dir(3,2)); L_tool];
    P4_m = [r2*cos(limb_dir(4,2)); r2*sin(limb_dir(4,2)); L_tool+h];
    P5_m = [r2*cos(limb_dir(5,2)); r2*sin(limb_dir(5,2)); L_tool+h];
    P_m = [P1_m P2_m P3_m P4_m P5_m];

    basic_paras.l_max = l_max;
    basic_paras.l_min = l_min;
    basic_paras.R1 = R1;
    basic_paras.R2 = R2;
    basic_paras.H = H;
    basic_paras.r1 = r1;
    basic_paras.r2 = r2;
    basic_paras.h = h;
    basic_paras.L_tool = L_tool;
    basic_paras.limb_dir = limb_dir;
    basic_paras.B = B;
    basic_paras.P_m = P_m;
    basic_paras.l0_seq = l0_seq;
    basic_paras.joint_u_angle_tilt = joint_u_angle_tilt;
    basic_paras.unit_para = unit_para;

end