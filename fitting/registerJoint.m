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

function template = registerJoint(scan,template,saveDir,modelDir,initDir)
fprintf('registerJoint()\n');

if (nargin < 5)
    initDir = saveDir;
end

templateNew = [];

% fit shape and pose parameters to each scan independently
for i=1:length(scan)
    % rigid alignment of template to the scan
    scan(i).T = [];
    scan(i) = rigidAlignTemplate2Scan(scan(i),template(i));
    visFit(scan(i),template(i));
    visLandmarks(scan(i),template(i));
    tmpl = fitModel(scan(i),template(i),modelDir,initDir{i});
    if (isempty(templateNew))
        templateNew = tmpl;
    else
        templateNew(i) = tmpl;
    end
end
template = templateNew;

shapeParams = zeros(1,length(template(1).shapeParams));

% average shape parameters, use independent pose parameters
for i=1:length(template)
    shapeParams = shapeParams + template(i).shapeParams;
end

for i=1:length(template)
    template(i).shapeParams = shapeParams ./length(template);
end

% fit shape and pose parameters
try
    for i=1:length(template)
        fname = [saveDir{i} '/poseShapeJoint'];
        fprintf('load %s\n',fname);
        load(fname,'points');
        template(i).points = points;
    end
catch
    tic
    template = fitPoseShapeJoint(scan,template,modelDir);
    for i=1:length(template)
        fname = [saveDir{i} '/poseShapeJoint'];
        points4     = [template(i).points ones(size(template(i).points,1),1)]/scan(i).T';
        points      = points4(:,1:3);
        poseParams  = template(i).poseParams;
        shapeParams = template(i).shapeParams;
        pointsIdxsScanNN = template(i).pointsIdxsScanNN;
        pointsHasValidNN = template(i).pointsHasValidNN;
        T = scan(i).T;
        dist = template(i).dist;
        save(fname,'points','poseParams','shapeParams','pointsIdxsScanNN','pointsHasValidNN','dist','T');
    end
    toc
end

end