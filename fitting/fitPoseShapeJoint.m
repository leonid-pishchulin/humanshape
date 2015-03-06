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

function template = fitPoseShapeJoint(scan,template,modelDir,threshNormAngle)

fprintf('fitPoseShapeJoint()\n');

if (nargin < 4)
    threshNormAngle = 60; % degrees
end;

% load idxHand
load('VertexIdxSpecParts', 'idxHand');
load([modelDir '/evectors'], 'evectors');
assert(template(1).nPCA <= size(evectors,1));
evectors = evectors(1:template(1).nPCA,:);

err = 0;
isValidNN = zeros(length(template(1).points),length(scan));

assert(length(scan) == length(template));

% loop over scans
for i = 1:length(scan)
    
    poseParams = template(i).poseParams;
    
    % optimize over global scale during pose and shape fitting
    sc = poseParams(end);
    shapeParams = [template(i).shapeParams sc];
    
    pointsSM = changeShapePose(poseParams(1:end-1), shapeParams(1:end-1),evectors,modelDir);
    pointsSM = pointsSM * sc;
    
    %% find normals shape model
    nPoints = size(pointsSM,1);
    normalsSM = getNormals(pointsSM,template(i).faces);
    normalsSM = getNormals1Face(1:nPoints,template(i).faces,normalsSM);
    
    %% compute NN
    [idxsNN, distNN] = knnsearch(scan(i).points,pointsSM);
    
    template(i).pointsIdxsScanNN = idxsNN;
    
    %% find normals scan
    assert(size(idxsNN,1) == nPoints);
    normalsScanAll = getNormals1Face(scan(i).pointsIdxs,scan(i).faces,scan(i).normals);
    normalsScan = normalsScanAll(idxsNN,:);
    
    %% check the angle between normals
    isValidNN(:,i) = checkAngle(normalsScan,normalsSM,threshNormAngle);
    
    % do not register open to closed hands
    isValidNN(idxHand,i) = 0;
    
    if (~isempty(template(i).idxsUse))
        % use only subset of vertices
        isValidNN(:,i) = isValidNN(:,i).*template(i).idxsUse;
    end
    
    distAll = distNN' * isValidNN(:,i);
    err = err + distAll/sum(isValidNN(:,i));
end

err = err ./ length(scan);

eps_err = 1e-3;
errPrev = err+eps_err+1;

%% optimization loop
while abs(errPrev - err) > eps_err    
    
    errPrev = err;
    
    fprintf('fit pose\n');
    
    % fit pose to each scan separately
        
    shapeParams = [template(1).shapeParams template(1).poseParams(end)];
    scAll = ones(length(scan),1);
    for i = 1:length(scan)
        fprintf('scan: %d\n',i);
        assert(sum(template(i).shapeParams - shapeParams(1:end-1))==0);
        shapeParams(end) = template(i).poseParams(end);
        % optimization parameters
        [options,poseLB,poseUB,shapeLB,shapeUB] = getOptionsOptimizer(modelDir,template(i));
        
        poseLB(template(i).poseParamsIgnoreIdxs) = 0;
        poseUB(template(i).poseParamsIgnoreIdxs) = 0;
        
        poseParams = template(i).poseParams;
        pointsScan = scan(i).points;
        isValidNNscan = isValidNN(:,i);
        idxsNN = template(i).pointsIdxsScanNN;
        
%         visFit(scan(i),template(i));
        
        [poseParams, ~] = fmincon(@PoseFunc,poseParams,[],[],[],[],poseLB,poseUB,[],options);
        [pointsSM, ~] = changeShapePose(poseParams(1:end-1),shapeParams(1:end-1),evectors,modelDir);
        sc = poseParams(end);
        pointsSM = sc * pointsSM;
        idxsNN = knnsearch(scan(i).points,pointsSM);
        
        % check the angle between normals
        normalsScanAll = getNormals1Face(scan(i).pointsIdxs,scan(i).faces,scan(i).normals);
        normalsScan = normalsScanAll(idxsNN,:);
        normalsSM = getNormals(pointsSM, template(i).faces);
        normalsSM = getNormals1Face(1:nPoints,template(i).faces,normalsSM);
        isValidNN(:,i) = checkAngle(normalsScan,normalsSM,threshNormAngle);
        
        % do not register open to closed hands
        isValidNN(idxHand,i) = 0;
        
        if (~isempty(template(i).idxsUse))
            % use only subset of vertices
            isValidNN(:,i) = isValidNN(:,i).*template(i).idxsUse;
        end
        template(i).points = pointsSM;
        template(i).poseParams = poseParams;
        template(i).pointsIdxsScanNN = idxsNN;
        scAll(i) = sc;
        
        % check
        d = sqrt(sum((scan(i).points(idxsNN, :) - template(i).points) .^ 2, 2));
        E(i) = sum(d.*isValidNN(:,i));
    end
    
    fprintf('fit shape\n');
    
    sc = mean(scAll);
    shapeParams(end) = sc;
    err = 0;
    
    [~,~,~,shapeLB1,shapeUB1] = getOptionsOptimizer(modelDir,template(1));
    [options,~,~,shapeLB2,shapeUB2] = getOptionsOptimizer(modelDir,template(2));
    assert(sum(shapeLB1 ~= shapeLB2) == 0);
    assert(sum(shapeUB1 ~= shapeUB2) == 0);
    
    [shapeParams, ~] = fmincon(@ShapeFunc,shapeParams,[],[],[],[],shapeLB1,shapeUB1,[],options);
    sc = shapeParams(end);
    for i = 1:length(template)
        poseParams = [template(i).poseParams shapeParams(end)];
        [pointsSM, ~] = changeShapePose(poseParams(1:end-1),shapeParams(1:end-1),evectors,modelDir);
        % update shape parameters (the same for both templates)
        template(i).shapeParams = shapeParams(1:end-1);
        % update global scaling factor
        template(i).poseParams(end) = sc;
        pointsSM = sc * pointsSM;
        [idxsNN, distNN] = knnsearch(scan(i).points,pointsSM);
        template(i).points = pointsSM;
        
        % check the angle between normals
        normalsScanAll = getNormals1Face(scan(i).pointsIdxs,scan(i).faces,scan(i).normals);
        normalsScan = normalsScanAll(idxsNN,:);
        normalsSM = getNormals(pointsSM,template(i).faces);
        normalsSM = getNormals1Face(1:nPoints,template(i).faces,normalsSM);
        isValidNN(:,i) = checkAngle(normalsScan,normalsSM,threshNormAngle);
        
        % do not register open to closed hands
        isValidNN(idxHand,i) = 0;
        
        if (~isempty(template(i).idxsUse))
            % use only subset of vertices
            isValidNN(:,i) = isValidNN(:,i).*template(i).idxsUse;
        end
        
        template(i).pointsIdxsScanNN = idxsNN;
        template(i).pointsHasValidNN = isValidNN(:,i);
        
        distAll = distNN' * isValidNN(:,i);
        err = err + distAll/sum(isValidNN(:,i));
        
        template(i).dist = distNN;
    end
    err = err ./ length(template);
end

fprintf('done\n');

    %% fit pose objective function
    function E = PoseFunc(poseParam)
        scTmp = poseParam(1,32);
        [pointsSM, ~] = changeShapePose(poseParam(1:end-1),shapeParams(1:end-1),evectors,modelDir);
        pointsSM = scTmp * pointsSM;
        pointsScanNN = pointsScan(idxsNN, :);
        dist = sqrt(sum((pointsScanNN - pointsSM) .^ 2, 2));
        E = sum(dist.*isValidNNscan);
    end

    %% fit shape objective function
    function E = ShapeFunc(shapeParam)
        scTmp = shapeParam(end);
        E = 0;
        for n=1:length(scan)
            poseParams = template(n).poseParams;
            [pointsSM, ~] = changeShapePose(poseParams(1:end-1),shapeParam(1:end-1),evectors,modelDir);
            pointsSM = scTmp * pointsSM;
            idxsNN = template(n).pointsIdxsScanNN;
            pointsScanNN = scan(n).points(idxsNN, :);
            dist = sqrt(sum((pointsScanNN - pointsSM) .^ 2, 2));
            E = E + sum(dist.*isValidNN(:,n));
        end
        E = E/length(scan);
    end

end
