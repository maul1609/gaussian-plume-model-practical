% Overlay concentrations on map
load map_green_lane;
figure
imagesc(ximage,yimage,A)
set(gca,'ydir','normal') 
xlabel('x (m)');ylabel('y (m)');
h1=axes('position',get(gca,'position'));
[c,h]=contour(x,y,(mean(C1,3)).*1e6);
set(gca,'visible','off')
set(h,'linewidth',5);
clabel(c,h,'fontsize',20)
