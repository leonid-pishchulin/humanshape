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

function scan = rigidAlignTemplate2Scan(scan,template)

regParams = absor(template.points(template.landmarksIdxs(scan.landmarkIdxs),:)',scan.landmarks(scan.landmarkIdxs,:)','doScale',true);

% template.points = points';
% template.points = template.points*regParams.R*regParams.s;
T = eye(4);
T(1:3,1:3) = regParams.R/regParams.s;
pointsOrig = scan.points;

points = [scan.points ones(size(scan.points,1),1)]*T;
scan.points = points(:,1:3);
t = mean(template.points) - mean(scan.points);
scan.points = scan.points + repmat(t,size(scan.points,1),1);

T(4,1:3) = t;
T = T';

landmarks = [scan.landmarks, ones(size(scan.landmarks,1),1)]*T';
scan.landmarks = landmarks(:,1:3);

points = [pointsOrig ones(size(scan.points,1),1)]*T';
scan.points = points(:,1:3);

% update 16.01.2015
if (isfield(scan,'normals'))
    Tnew = T;
    Tnew(1:3,4) = 0;
    normals = [scan.normals ones(size(scan.normals,1),1)]*Tnew';
    scan.normals = normals(:,1:3);
end

scan.T = T;
% [alpha,theta,gamma] = rotationMatrix(regParams.R);
% template.poseParams(1:3) = [-alpha+pi,-pi+theta,-gamma+pi];
% template.poseParams(1:3) = [0 0 0];
% template.poseParams(4:6) = t;

end