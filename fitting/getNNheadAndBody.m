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

function [pointsIdxsScanNNall] = getNNheadAndBody(scan,template)
% scan Nearest Neighbour points for template points
% treat head and body separately

% distance threshold for NN
% NNdistThresh = 20;
% for bootstrap > 2
NNdistThresh = 50;

NNdistThreshHead = 50;

% maximum difference between normals of NN for the head
maxHeadNormalAngle = 30;

% maximum difference between normals of NN for the body
maxBodyNormalAngle = 60;

nPointsTemplate = size(template.points,1);

[templateIdxAbove, templateIdxBelow] = getHeadTopIdx(template.points,template.points);
scanIdxAbove                         = getHeadTopIdx(template.points,scan.points);
scanIdxAbove = intersect(scan.pointsIdxsFaces,scanIdxAbove);

% Finding the appropriate Nearest Neighbour for the above cutting plane points
pointsTemplateAbove = template.points(templateIdxAbove,:);
pointsScanAbove     = scan.points    (scanIdxAbove,:);

nPointsTemplateAbove = length(templateIdxAbove);
nPointsScanAbove     = length(scanIdxAbove);

templateNormalsAbove = template.normals(templateIdxAbove,:);

pointsIdxsScanNNall = nan(nPointsTemplate,1);

for i = 1:nPointsTemplateAbove
    diff = pointsScanAbove - repmat(pointsTemplateAbove(i,:),[nPointsScanAbove 1]);
    normal = repmat(templateNormalsAbove(i,:),[nPointsScanAbove 1]);
    
    isAligned = checkAngle(diff,normal,maxHeadNormalAngle);
    
    dist = sqrt(sum(diff.^2,2));
    idxs = find(isAligned > 0 & dist <= NNdistThreshHead);
    
    if (~isempty(idxs))
        [val, idx] = min(dist(idxs));
        pointsIdxsScanNNall(templateIdxAbove(i)) = scanIdxAbove(idxs(idx));
    end
end

[pointsIdxsScanNN, NNdist] = knnsearch(scan.points,template.points);

normalsScanNN = scan.normals1face(pointsIdxsScanNN,:);
isAligned = checkAngle(normalsScanNN,template.normals,maxBodyNormalAngle);

idxs = intersect(find(isAligned > 0), templateIdxBelow);

% set distance threshold
idxs2 = find(NNdist <= NNdistThresh);
idxs = intersect(idxs, idxs2);

pointsIdxsScanNNall(idxs) = pointsIdxsScanNN(idxs);

end