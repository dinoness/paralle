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

%初始位置时 动平台相对静平台的位姿
Q=eye(3);
%Q=rotz(0)*roty(0)*rotx(45);

%支链球体移动后的球心坐标
sphere=[xa;ya;za]-Q*[xb;yb;zb];
u=sphere(1,:);
v=sphere(2,:);
w=sphere(3,:);

zmin=min(w);%Z方向工作空间最小高度
zmax=max(w)+max(roumax);%Z方向工作空间最大高度

discreteNum=10;%每条弧的离散点数量
numCir=size(sphere,2);%支链个数

% Rot=VectorRotMartix([1,0,0],[0,1,0]) 切割面法向量与XOY平面法向量 的变换矩阵
tic
%XY平面切割遍历
for range=zmin:4:zmax
    %平移矩阵
    Tra=transl(0,0,-range);
    %旋转矩阵
    Rot=eye(4);
    %齐次变换矩阵
    P=Tra*Rot;
    %构建球的齐次坐标 并变换到XOY平面
    sphere_P=P*[sphere;ones(1,numCir)];
    %求平面内的重合弧
    arcList=TraversePlane(sphere_P,[roumin;roumax]);
    numArc=size(arcList,2);%这个平面中 重合弧的数量
    if numArc==0
        continue;
    end
    %离散每条重合弧
    points=ones(4,numArc*discreteNum);%记录该平面内的离散点 齐次坐标
    for i = 1:numArc
        xr = arcList(1,i);
        yr = arcList(2,i);
        r = arcList(3,i);
        startAngle = arcList(4,i);
        endAngle = arcList(5,i)+startAngle;
        theta = linspace(startAngle,endAngle,discreteNum);
        x = xr+r*cosd(theta);
        y = yr+r*sind(theta);
        %XY平面内的离散点 齐次坐标
        pointsXY(1:2,:)=[x;y];
        pointsXY(3,:)=0;
        pointsXY(4,:)=1;
        %XY平面中的点反变换回原切割平面
        pointsInitial=inv(P)*pointsXY;
        plot3(pointsInitial(1,:),pointsInitial(2,:),pointsInitial(3,:));
        hold on;
    end
end
toc

tic
%绕Z轴旋转切割遍历
for range=0:10:360
    %平移矩阵
    Tra=transl(0,0,0);
    %旋转矩阵
    OA=[sind(range),cosd(range),0];
    OB=[0,0,1];
    Rot=VectorRotMartix(OA,OB);%OA为切割面法向量 OB为XOY平面法向量
    %齐次变换矩阵
    P=Tra*Rot;
    %构建球的齐次坐标 并变换到XOY平面
    sphere_P=P*[sphere;ones(1,numCir)];
    %求平面内的重合弧
    arcList=TraversePlane(sphere_P,[roumin;roumax]);
    numArc=size(arcList,2);%这个平面中 重合弧的数量
    if numArc==0
        continue;
    end
    %离散每条重合弧
    points=ones(4,numArc*discreteNum);%记录该平面内的离散点 齐次坐标
    for i = 1:numArc
        xr = arcList(1,i);
        yr = arcList(2,i);
        r = arcList(3,i);
        startAngle = arcList(4,i);
        endAngle = arcList(5,i)+startAngle;
        theta = linspace(startAngle,endAngle,discreteNum);
        x = xr+r*cosd(theta);
        y = yr+r*sind(theta);
        %XY平面内的离散点 齐次坐标
        pointsXY(1:2,:)=[x;y];
        pointsXY(3,:)=0;
        pointsXY(4,:)=1;
        %XY平面的点反变换回原始平面
        pointsInitial=inv(P)*pointsXY;
        %找到工作空间中Z轴于0的列 的索引
        LessZero=find(pointsInitial(3,:)<0);
        %删除Z轴小于0的数据 对于工作空间没有意义
        pointsInitial(:,LessZero)=[];
        plot3(pointsInitial(1,:),pointsInitial(2,:),pointsInitial(3,:));
        hold on;
    end
end
toc
view([6,2,2]);
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
