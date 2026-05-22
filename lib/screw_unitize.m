function x = screw_unitize(x, screw_type)
%SCREW_UNITIZE 对旋量/力旋量做统一归一化，便于显示和后续比较
% twist  : 优先按角速度部分归一；若角速度接近0，则按线速度部分归一
% wrench : 优先按力部分归一；若力接近0，则按力矩部分归一

    tol = 1e-12;

    switch lower(screw_type)
        case 'twist'
            if norm(x(1:3)) > tol
                x = x / norm(x(1:3));
            elseif norm(x(4:6)) > tol
                x = x / norm(x(4:6));
            end

        case 'wrench'
            if norm(x(1:3)) > tol
                x = x / norm(x(1:3));
            elseif norm(x(4:6)) > tol
                x = x / norm(x(4:6));
            end

        otherwise
            error('screw_unitize: 未知 screw_type = %s', screw_type);
    end
end
