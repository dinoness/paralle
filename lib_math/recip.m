function x = recip(s1,s2)
    % 计算两个旋量的互易积
s11 = s1(1:3);
s12 = s1(4:6);
s21 = s2(1:3);
s22 = s2(4:6);
x = dot(s11,s22) + dot(s12,s21);
end