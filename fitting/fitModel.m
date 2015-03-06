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

function template = fitModel(scan,template,modelDir,initDir)

% fit pose parameters
if (~isempty(scan.landmarks))
    fname = [initDir '/poseInit'];
    try
        fprintf('load %s\n',fname);
        load(fname,'poseParams');
        template.poseParams = poseParams;
    catch
        tic
        [template,dist] = fitPose(scan,template,modelDir);
        % visualize fitted landmarks
        visLandmarks(scan,template);
        toc
        % transform points into the scan's coordinate system
        poseParams = template.poseParams;
        points4    = [template.points ones(size(template.points,1),1)]/scan.T';
        points     = points4(:,1:3);
        T = scan.T;
        save(fname,'points','poseParams','dist','T');
    end
end

% fit pose and shape parameters
fname = [initDir '/poseShapeInit'];
try
    fprintf('load %s\n',fname);
    load(fname,'points','poseParams','shapeParams','pointsIdxsScanNN','pointsHasValidNN');
    points4     = [points ones(size(template.points,1),1)]*scan.T';
    points      = points4(:,1:3);
    template.points = points;
    template.poseParams = poseParams;
    template.shapeParams = shapeParams;
    template.pointsIdxsScanNN = pointsIdxsScanNN;
    template.pointsHasValidNN = pointsHasValidNN;
catch
    tic
    [template,dist] = fitPoseShape(scan,template,modelDir);
    % visualize fitted mesh
    visFit(scan,template);
    % transform points into the scan's coordinate system
    points4     = [template.points ones(size(template.points,1),1)]/scan.T';
    points      = points4(:,1:3);
    poseParams  = template.poseParams;
    shapeParams = template.shapeParams;
    pointsIdxsScanNN = template.pointsIdxsScanNN;
    pointsHasValidNN = template.pointsHasValidNN;
    T = scan.T;
    save(fname,'points','poseParams','shapeParams','pointsIdxsScanNN','pointsHasValidNN','dist','T');
    toc
end
end