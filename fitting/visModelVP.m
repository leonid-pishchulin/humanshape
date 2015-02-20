function visModelVP(expidx,VP,idxShape,sign,bSave)

p = expParams(expidx);

modelDir = p.modelDir;
if (bSave)
    figuresDir = p.figuresDir;
end

load([p.modelInDir '/evalues'], 'evalues');
load([p.modelInDir '/evectors'], 'evectors');
load('facesShapeModel.mat','faces');
shapeParams = zeros(1,20);

shapeParams(idxShape) = sign*3*sqrt(evalues(idxShape));
points = changeShapePose(zeros(1,31), shapeParams, evectors, p.modelInDir);
points(:,3) = points(:,3) - min(points(:,3));

template.points = points;
template.faces = faces;

fontSize = 18;

figure(100); clf;
set(0,'DefaultAxesFontSize', fontSize)
set(0,'DefaultTextFontSize', fontSize)

Xl = [-250 250];
Yl = [-250 250];
Zl = [0 2100];

showmodel(template.points,template.faces,p.colorName,[],0);hold on;
set(gca,'XLim',Xl,'YLim',Yl,'ZLim',Zl,'FontSize',fontSize);

view(VP,0);
axis equal;
axis off;
grid on;

[path,name] = fileparts(p.modelInDir);
[~,name2] = fileparts(path);

if (bSave)
    print(gcf,'-dpng',[figuresDir '/' name2 '-idxShape-' num2str(idxShape) '-' num2str(sign) '-vp-' num2str(VP) '.png']);
    printpdf([figuresDir '/' name2 '-idxShape-' num2str(idxShape) '-' num2str(sign) '-vp-' num2str(VP) '.pdf']);
end

end