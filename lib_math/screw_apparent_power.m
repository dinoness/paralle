function [P_app, P_real, zeta, info] = screw_apparent_power(S1, S2, varargin)
% 计算 S1 与 S2 之间的视在功率
% =======================约定：S1是力，S2是位移=======================
%
% 输入:
%   S1 : 6x1 运动旋量 [omega; v]
%   S2 : 6x1 力旋量   [f; m]
%
% 可选输入:
%   'dMax', value : 指定文献公式中的 d_max
%                   若不指定，则默认用当前两旋量轴线之间的公垂距 d
%   'Tol', value  : 判断零向量的阈值，默认 1e-12
%
% 输出:
%   P_app  : 视在功率
%   P_real : 实际功率，即互易积绝对值
%   zeta   : 能效系数 P_real / P_app
%   info   : 中间变量，包括 h1, h2, d, dMax, kind1, kind2

tol = 1e-12;
dMax_user = [];
p_cstr = [];

for k = 1:2:numel(varargin)
    switch lower(varargin{k})
        case 'dmax'
            dMax_user = varargin{k+1};
        case 'tol'
            tol = varargin{k+1};
        case 'p_cstr'
            p_cstr = varargin{k+1};
        otherwise
            error('Unknown option: %s', varargin{k});
    end
end

S1 = S1(:);
S2 = S2(:);

if numel(S1) ~= 6 || numel(S2) ~= 6
    error('S1 and S2 must both be 6x1 screw vectors.');
end

Omega = [zeros(3,3), eye(3);
         eye(3), zeros(3,3)];

[S1u, h1, q1, kind1] = normalize_screw(S1, tol);
[S2u, h2, q2, kind2] = normalize_screw(S2, tol);

P_real_signed = S2u.' * Omega * S1u;
P_real = abs(P_real_signed);

if strcmp(kind1, 'finite') && strcmp(kind2, 'finite')
    % d = common_normal_distance(S1u(1:3), q1, S2u(1:3), q2, tol);
    if isempty(p_cstr)
        error('p_cstr is empty while calculate d_max.');
    else
        d = potential_max_distance(S1u(1:3), p_cstr, S2u(1:3), q2, tol);
    end

    if isempty(dMax_user)
        dMax = d;
    else
        dMax = dMax_user;
    end

    P_app = sqrt((h1 + h2)^2 + dMax^2);

elseif strcmp(kind1, 'infinite') && strcmp(kind2, 'finite')
    v1 = S1u(4:6);
    f2 = S2u(1:3);

    P_app = norm(v1) * norm(f2);
    d = NaN;
    dMax = NaN;

elseif strcmp(kind1, 'finite') && strcmp(kind2, 'infinite')
    omega1 = S1u(1:3);
    m2 = S2u(4:6);

    P_app = norm(omega1) * norm(m2);
    d = NaN;
    dMax = NaN;

else
    P_app = 0;
    d = NaN;
    dMax = NaN;
end

if P_app > tol
    zeta = P_real / P_app;
else
    zeta = 0;
end

info.S1_unit = S1u;
info.S2_unit = S2u;
info.h1 = h1;
info.h2 = h2;
info.q1 = q1;
info.q2 = q2;
info.d = d;
info.dMax = dMax;
info.kind1 = kind1;
info.kind2 = kind2;
info.P_real_signed = P_real_signed;

end

function [Su, h, q, kind] = normalize_screw(S, tol)

a = S(1:3);
b = S(4:6);

if norm(a) > tol
    Su = S / norm(a);  % 依据原部矢量进行单位化
    a = Su(1:3);
    b = Su(4:6);

    h = dot(a, b);  % 是否要除分母？除不除无所谓
    q = cross(a, b);
    kind = 'finite';

elseif norm(b) > tol
    Su = S / norm(b);
    h = Inf;
    q = [NaN; NaN; NaN];
    kind = 'infinite';

else
    error('Invalid screw: both parts are near zero.');
end

end

function d = common_normal_distance(a1, q1, a2, q2, tol)

a1 = a1 / norm(a1);
a2 = a2 / norm(a2);

n = cross(a1, a2);

if norm(n) > tol
    d = abs(dot(q2 - q1, n)) / norm(n);
else
    d = norm(cross(q2 - q1, a1));
end


end

function d = potential_max_distance(a1,p_constrain,a2,q2,tol)
    l = q2 - p_constrain;
    d = norm(cross(a2, l));
end

