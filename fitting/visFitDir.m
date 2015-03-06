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

function res = visFitDir(fitDir,fitType,bTransform,angle)

if (nargin < 2)
    fitType = 2; % NRD fit
end
if (nargin < 3)
    bTransform = false; % transform the point into SM coordinate system
end
if (nargin < 4)
    angle = [-45 45]; % viewpoint angle
end

try
    switch fitType
        case 0
            load([fitDir '/poseInit'],'points','T');
        case 1
            load([fitDir '/poseShapeInit'],'points','T');
        case 2
            load([fitDir '/NRD'],'points','pointsIdxsScanNN');
        case 3
            load([fitDir '/poseShapeJoint.mat'],'points','T');
    end
catch
    res = 0;
    return;
end

template.points = points;
load([fitDir '/scan'],'points');
scan.points = points;

if (bTransform)
    points4     = [template.points ones(size(template.points,1),1)]*T';
    template.points  = points4(:,1:3);
    points4     = [scan.points ones(size(scan.points,1),1)]*T';
    scan.points  = points4(:,1:3);
end

load('facesShapeModel','faces');
template.faces = faces;

visFit(scan,template,angle);

res = 1;
end