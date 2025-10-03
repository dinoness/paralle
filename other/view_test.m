% clear; close all;

% t = 0:pi/10:10*pi;
% x = sin(t);
% y = cos(t);

% for ii = 2:length(t)
%     figure(1); 
%     % plot points
%     p = scatter3(x(ii),y(ii),t(ii),'red','LineWidth',2); hold on;
    
%     % plot traces
%     plot3(x(ii-1:ii),y(ii-1:ii),t(ii-1:ii),'r-',"LineWidth",2);
    
%     % reverse Z axis
%     set(gca,'Zdir','reverse');
    
%     % set fixed figure size
%     set(gcf,'Position',[400 80 600 600]);

%     % set fixed XYZ axis
%     xlim([min(x),max(x)]);
%     ylim([min(y),max(y)]);
%     zlim([min(t),max(t)]);

%     % set fixed XYZ ticks
%     xticks([min(x),(max(x)+min(x))/2,max(x)]);
%     yticks([min(y),(max(y)+min(y))/2,max(y)]);
%     zticks([min(t),(max(t)+min(t))/2,max(t)]);
%     zticklabels({'0','5\pi','10\pi'});

%     % Freeze the aspect ratio
%     axis vis3d;
    
%     % rotate the camera, save current figure as tif
%     view(140 + 2*ii,15);
%     pause(0.1);
%     % saveas(gcf,['animation3D\',num2str(ii),'.tif']);

%     % delete point, keep trace
%     delete(p);
% end

clear; clc;
x = 0:0.01:1;
n = 3;
y = x.^n;
plot(x,y,'LineWidth',3)
title(['y = x^n,  n = ' num2str(n) ])

n = 1:0.5:5;
nImages = length(n);
fig = figure;
for idx = 1:nImages
    y = x.^n(idx);
    plot(x,y,'LineWidth',3)
    title(['y = x^n,  n = ' num2str( n(idx)) ])
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);
end
close;

figure;
for idx = 1:nImages
    subplot(3,3,idx)
    plot(x, x.^idx,'LineWidth',3)
    title(['y = x^', num2str(idx) ])
end
filename = 'testAnimated.gif'; % Specify the output file name
for idx = 1:nImages
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.3);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
    end
end