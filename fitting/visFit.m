function visFit(scan,template,angle)

if (nargin < 3)
    angle = [-45 45];
end

clf;
offset = 25;
Xl = [min(template.points(:,1))-offset max(template.points(:,1))+offset];
Yl = [min(template.points(:,2))-offset max(template.points(:,2))+offset];
Zl = [min(template.points(:,3))-offset max(template.points(:,3))+offset];

if (length(angle) == 2)
    subplot(1,2,1);hold on;
    showmodel(template.points,template.faces,'r',[],0);
    plot3(scan.points(:,1),scan.points(:,2),scan.points(:,3),'b.','MarkerSize',2);
    set(gca,'XLim',Xl,'YLim',Yl,'ZLim',Zl);
    view(angle(1),0);
    axis equal;
    
    subplot(1,2,2);hold on;
    showmodel(template.points,template.faces,'r',[],0);
    plot3(scan.points(:,1),scan.points(:,2),scan.points(:,3),'b.','MarkerSize',2);
    set(gca,'XLim',Xl,'YLim',Yl,'ZLim',Zl);
    view(angle(2),0);
    axis equal;
else
    showmodel(template.points,template.faces,'r',[],0);
    hold on;
    showmodel(scan.points,template.faces,'b',[],0);
    set(gca,'XLim',Xl,'YLim',Yl,'ZLim',Zl);
    view(angle,0);
    axis equal;
end
%     fname = [fitDir '/' name{i} '.png'];
%     print(gcf,'-dpng',fname);
  % close(100);
end