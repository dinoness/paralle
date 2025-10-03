clear all

%静平台 坐标参数
xa=[92.58 132.58 40 -40 -132.58 -92.58];
ya=[99.64 30.36 -130 -130 30.36 99.64];
za=23.1*ones(1,6);
%动平台 坐标参数
xb=[30 78.22 48.22 -48.22 -78.22 -30];
yb=[73 -10.52 -62.48 -62.48 -10.52 73];
zb=-37.1*ones(1,6);
%支链长度范围
roumin=454.5*ones(1,6);
roumax=504.5*ones(1,6);


% XOY切割
figure(1);

seq_z = 460 : 20 : 520;
for iz = 1 : length(seq_z)
    z_slice = seq_z(iz);
    %内外圆
    circle_in = [xa;ya;sqrt(roumin.^2 - (za - z_slice).^2);ones(1,length(xa)).*(-1)]; %内圆
    circle_out = [xa;ya;sqrt(roumax.^2 - (za - z_slice).^2);ones(1,length(xa))]; %外圆

    circles = [circle_in circle_out]; %所有圆
    ArcDraw(circle_out, z_slice);

end
