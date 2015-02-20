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