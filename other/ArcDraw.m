function ArcDraw(circle, z)

    [~,col] = size(circle);
    %-------交点列表 用于画图-------
    crossPointsList = [];
    %---------弧列表初始化----------
    arcList = [];  % 弧列表 x y r startAngle arcAngle
    %---------圆列表初始化----------
    cirList(col) = Circle();%建立类对象数组
    for i = 1:col
        cirList(i).cirId = i;
        cirList(i).cirPro = circle(1:3,i) ;
        cirList(i).cirType = circle(4,i) ;
    end
    %---------遍历所有圆，两两求交--------
    for i = 1:col
        for j = i+1:col
            %两圆交点
            crossPoint = ArcIntersection(cirList(i),cirList(j));
            %如果有交点
            if ~isempty(crossPoint)
                %将交点放入交点列表
                crossPointsList(:,end+1:end+2) = crossPoint;
            end
        end
        PlotCir(cirList(i), z);%画圆
        hold on;
        cirList(i).CrossPointsAngle();
        %画出两交点的中心点
        if ~isempty(cirList(i).middlePoints)
            % plot(cirList(i).middlePoints(1,:), ...
            %     cirList(i).middlePoints(2,:), ...
            %     'b.','MarkerSize',25);
            hold on;
        end
    end
    %--------画出交点-----------
    if isempty(crossPointsList)
        fprintf('无任何交点\n');
    else
        % plot(crossPointsList(1,:),crossPointsList(2,:),'r.','MarkerSize',25);
        hold on;
    end

    % 判断哪些弧属于重合弧
    for i = 1:col
        arcFlag = (1:cirList(i).cirArcNum);
        for j = 1:col
            %跳过自身判断
            if (i == j)
                continue
            end
            pos = cirList(i).middlePoints - cirList(j).cirPro(1:2,:);
            d = (pos(1,:).^2+pos(2,:).^2).^0.5;
            %如果是外圆
            if cirList(j).cirType == 1 
                arcOder_temp = find(d <= cirList(j).cirPro(3)); %满足要求的弧的序号
            else %如果是内圆
                arcOder_temp = find(d >= cirList(j).cirPro(3)); %满足要求的弧的序号
            end
            arcFlag = intersect(arcFlag,arcOder_temp);
        end
        %如果没有n次相交弧
        if isempty(arcFlag)
            fprintf('%d号圆没有n次相交弧\n', cirList(i).cirId);
            continue;
        else
            arcList_temp = ones(3,size(arcFlag,2)).*cirList(i).cirPro;
            arcList_temp(4,:) = cirList(i).cirArcList(1,arcFlag);
            arcList_temp(5,:) = cirList(i).cirArcList(2,arcFlag);
            arcList = [arcList arcList_temp];
        end
    end
    %画圆弧
    % size(arcList)
    PlotArc(arcList, z)
    hold on
end


% 求两圆交点
function output = ArcIntersection(circ1, circ2)
    x1 = circ1.cirPro(1);
    y1 = circ1.cirPro(2);
    r1 = circ1.cirPro(3);
    x2 = circ2.cirPro(1);
    y2 = circ2.cirPro(2);
    r2 = circ2.cirPro(3);
    dis = ((x2-x1)^2+(y2-y1)^2)^0.5;

    if dis >= r1 + r2 % 两圆相离 无交点
        output = [];
        return;
    elseif dis <= abs(r1 - r2) % 内含 无交点
        output = [];
        return;
    end
    % 两圆心连线与X轴夹角，逆时针为正
    theta = atan2d(y2-y1,x2-x1);
    ad = (r1^2-r2^2+dis^2)/(2*dis);
    % 交点与圆心1连线的向量与两圆心连线的夹角
    alpha = acosd(ad/r1);
    crosspoint1 = [x1+r1*cosd(theta+alpha);y1+r1*sind(theta+alpha)];
    crosspoint2 = [x1+r1*cosd(theta-alpha);y1+r1*sind(theta-alpha)];

    circ1.cirCrossPoint(:,end+1:end+2) = [crosspoint1 crosspoint2];
    circ2.cirCrossPoint(:,end+1:end+2) = [crosspoint1 crosspoint2];
    output = [crosspoint1 crosspoint2];
end

% 画圆
function PlotCir(circ, z)
    theta = linspace(0,2*pi,100);
    x = circ.cirPro(3)*cos(theta)+circ.cirPro(1);
    y = circ.cirPro(3)*sin(theta)+circ.cirPro(2);
    plot3(x,y,ones(length(x))*z,'k-','LineWidth',0.5);
end

% 画圆弧
function PlotArc(arcList, z)
    for i = 1 : length(arcList(1,:))
        startAngle = arcList(4,i);
        arcAngle = arcList(5,i);
        angle = linspace(deg2rad(startAngle), deg2rad(startAngle+arcAngle), 50);
        x = arcList(3,i)*cos(angle)+arcList(1,i);
        y = arcList(3,i)*sin(angle)+arcList(2,i);
        plot3(x,y,ones(length(x))*z,'-','LineWidth',2,'Color',"#6495ED");
        hold on;
    end
end