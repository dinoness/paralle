function T = pos2trans(pos, B)
    % pos: x,y,z,phi(xoy平面内，与x轴夹角),theta(与z轴夹角)
    t = pos(1 : 3);
    
    phi = pos(4) / 180 * pi;
    theta = pos(5) / 180 * pi;

    z_axis = [sin(theta)*cos(phi);
         sin(theta)*sin(phi);
         cos(theta)];
    ObB1 = B(:, 1) - t;
    x_axis = cross(ObB1, z_axis) / norm(cross(ObB1, z_axis));
    y_axis = cross(z_axis, x_axis);
    R = [x_axis y_axis z_axis];

    T = [R t; 0 0 0 1];
end