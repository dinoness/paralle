classdef Circle < handle

    properties
        cirId  % 圆编号
        cirPro  % 圆属性 x y r
        cirType  % 圆类型 1外环 -1内环 0默认
        cirCrossPoint  % 交点列表，交点的x y
        cirArcList  % 圆弧列表，起点角度，旋转角度，逆时针
        cirArcNum  % 圆弧数量
        middlePoints  % 交点中间点 x y angle
    end

    methods
        % 构造函数
        function obj = Circle(circle, type)
            if nargin == 0
                obj.cirPro = zeros(3, 1);
                obj.cirType = 0;
            elseif nargin == 1
                obj.cirPro = circle;
                obj.cirType = 0;
            elseif nargin == 2
                obj.cirPro = circle;
                obj.cirType = type;
            end
        end

        % 计算交点与圆心连线 与 x轴夹角
        function CrossPointsAngle(obj)
            xc = obj.cirPro(1);
            yc = obj.cirPro(2);
            r = obj.cirPro(3);
            if ~isempty(obj.cirCrossPoint)
                vector = [obj.cirCrossPoint(1,:) - xc;
                          obj.cirCrossPoint(2,:) - yc];
                angles = atan2d(vector(2,:), vector(1,:));
                angleOrder = CrossPointsSort(angles);
                obj.cirArcList(1,:) = angleOrder;  % 每段圆弧的起始角度
                obj.cirArcList(2,:) = [diff(angleOrder), 360 - angleOrder(end) + angleOrder(1)];  % 每段圆弧的弧长
                obj.cirArcNum = size(obj.cirArcList, 2);
                % 生成探测点
                middleAngle = obj.cirArcList(1,:) + obj.cirArcList(2,:)/2;
                obj.middlePoints = [xc + r.*cosd(middleAngle);
                                    yc + r.*sind(middleAngle)];
            else
                obj.cirArcList(1,:) = 0;
                obj.cirArcList(2,:) = 360;
                obj.cirArcNum = 1;
                % 生成探测点
                middleAngle = obj.cirArcList(1,:) + obj.cirArcList(2,:)/2;
                obj.middlePoints = [xc + r.*cosd(middleAngle);
                                    yc + r.*sind(middleAngle)];
                fprintf('Circle %d has no cross point!\n', obj.cirId);
            end
        end


    end

end

%圆的交点排序 交点与X轴夹角从小到大排序 输入的是角度序列
function angleOrder = CrossPointsSort(angle)
    array = find(angle<0);
    %将小于0°的角度 转为大于0°
    angle(array)=angle(array)+360;
    %对圆的交点列表进行排序 因为是行排序 用sortrows时需要转置
    tempAngle = sortrows(angle',1);
    angleOrder = tempAngle';
end

