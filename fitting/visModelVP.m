%{
    This file is part of the evaluation of the 3D human shape model as described in the paper:

    Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobalt and Bernt Schiele
    Building Statistical Shape Spaces for 3D Human Modeling
    ArXiv, March 2015

    Please cite the paper if you are using this code in your work.
    
    Author: Leonid Pishchulin.

    The code may be used free of charge for non-commercial and
    educational purposes, the only requirement is that this text is
    preserved within the derivative work. For any other purpose you
    must contact the authors for permission. This code may not be
    redistributed without permission from the authors.
%}

function visModelVP(expidx,VP,idxShape,sign,bSave)

p = expParams(expidx);

if (bSave)
    figuresDir = p.figuresDir;
end

load([p.modelInDir '/evalues'], 'evalues');
load([p.modelInDir '/evectors'], 'evectors');
load('facesShapeModel.mat','faces');
evectors = evectors(1:p.nPCA,:);
evalues = evalues(1:p.nPCA);

shapeParams = zeros(1,20);
shapeParams(idxShape) = sign*3*sqrt(evalues(idxShape));
points = shapepose(zeros(1,31), shapeParams, evectors, p.modelInDir);
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