function Tad = adjoint_m(T)
    R = T(1:3, 1:3);
    t = T(1:3, 4);
    Tad = [R zeros(3,3); skew(t)*R R];
end