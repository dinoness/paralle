clear

R1 = 800;
r = 80;
Lmax = 1000;
Lmin = 650;
% L = 500;

L_seq = Lmin : 10 : Lmax;
% x_plant_seq = zeros(1, length(L_seq));
% h_plant_seq = zeros(1, length(L_seq));
x_plant_seq = [];
h_plant_seq = [];

for i = 1 : length(L_seq)
    L2 = L_seq(i);
    edge = sqrt(3) * R1;
    BD = Lmax * R1 / (R1 - r);
    AD = L2 * R1 / (R1 - r);
    DE = sqrt(BD^2 - (0.5*edge)^2);
    AE = 1.5 * R1;
    angle_AED = acos((AE^2 + DE^2 - AD^2) / (2 * AE * DE));
    h_D = DE * sin(angle_AED);
    h_plant = h_D * (R1 - r) / R1;
    x_D = 0 - 0.5 * R1 + DE * cos(angle_AED) ;
    x_plant = x_D * (R1 - r) / R1;
    
    x_plant_seq(end+1) = x_plant;
    h_plant_seq(end+1) = h_plant;
end
plot(x_plant_seq, h_plant_seq, 'LineWidth', 2);
axis equal
hold on

x_plant_seq = [];
h_plant_seq = [];
for i = 1 : length(L_seq)
        L1 = L_seq(i);
        edge = sqrt(3) * R1;
        BD = L1 * R1 / (R1 - r);
        AD = Lmax * R1 / (R1 - r);
        DE = sqrt(BD^2 - (0.5*edge)^2);
        AE = 1.5 * R1;
        angle_AED = acos((AE^2 + DE^2 - AD^2) / (2 * AE * DE));
        h_D = DE * sin(angle_AED);
        h_plant = h_D * (R1 - r) / R1;
        x_D = 0 - 0.5 * R1 + DE * cos(angle_AED) ;
        x_plant = x_D * (R1 - r) / R1;
        
        x_plant_seq(end+1) = x_plant;
        h_plant_seq(end+1) = h_plant;
end
plot(x_plant_seq, h_plant_seq, 'LineWidth', 2);





% edge = sqrt(3) * R1;
% BD = L * R1 / (R1 - r);
% AD = Lmax * R1 / (R1 - r);
% DE = sqrt(BD^2 - (0.5*edge)^2);
% AE = 1.5 * R1;
% angle_AED = acos((AE^2 + DE^2 - AD^2) / (2 * AE * DE));
% h_D = DE * sin(angle_AED);
% h_plant = h_D * (R1 - r) / R1;
% x_D = 0 - 0.5 * R1 + DE * cos(angle_AED) ;
% x_plant = x_D * (R1 - r) / R1;



