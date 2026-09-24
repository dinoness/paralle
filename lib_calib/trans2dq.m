function X = trans2dq(T)
%TRANS2DQ 齐次变换矩阵 → 单位对偶四元数 X = [qv; q0; pv; p0]
%   与 dq_fks_spr4ups 的约定一致：o = Im(p·q̄)，平移不取 1/2 因子。
%   旋转部分由 Shepperd 法从 R 提取，p = o ⊗ q（纯四元数左乘）。

    R = T(1:3, 1:3);
    o = T(1:3, 4);
    tr = trace(R);
    if tr > 0
        s = sqrt(tr + 1) * 2;
        q0 = s / 4;
        qv = [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)] / s;
    else
        [~, k] = max(diag(R));
        switch k
            case 1
                s = sqrt(1 + R(1,1) - R(2,2) - R(3,3)) * 2;
                q0 = (R(3,2)-R(2,3)) / s;
                qv = [s/4; (R(1,2)+R(2,1))/s; (R(1,3)+R(3,1))/s];
            case 2
                s = sqrt(1 + R(2,2) - R(1,1) - R(3,3)) * 2;
                q0 = (R(1,3)-R(3,1)) / s;
                qv = [(R(1,2)+R(2,1))/s; s/4; (R(2,3)+R(3,2))/s];
            case 3
                s = sqrt(1 + R(3,3) - R(1,1) - R(2,2)) * 2;
                q0 = (R(2,1)-R(1,2)) / s;
                qv = [(R(1,3)+R(3,1))/s; (R(2,3)+R(3,2))/s; s/4];
        end
    end
    if q0 < 0   % 统一符号，保证序列连续性
        q0 = -q0;  qv = -qv;
    end
    pv = q0*o + cross(o, qv);
    p0 = -dot(o, qv);
    X = [qv; q0; pv; p0];
end
