function [template,distNN] = fitPoseShape(scan,template,modelDir,threshNormAngle,bFitShape,bScale)
fprintf('\nfitPoseShape()\n');

if (nargin < 4)
    threshNormAngle = 60; % degrees
end;

if (nargin < 5)
    bFitShape = true; 
end;

if (nargin < 6)
    bScale = true;
end;

poseParams = template.poseParams;

% load idxHand
load('VertexIdxSpecParts', 'idxHand');

% optimize over global scale during pose and shape fitting
sc = poseParams(end);
shapeParams = [template.shapeParams sc];

load([modelDir '/evectors'], 'evectors');
assert(template.nPCA <= size(evectors,1));
evectors = evectors(1:template.nPCA,:);

pointsSM = shapepose(poseParams(1:end-1), shapeParams(1:end-1),evectors,modelDir);
pointsSM = pointsSM * sc;

%% find shape model mesh's normals
nPoints = size(pointsSM,1);
normalsSM = getNormals(pointsSM,template.faces);
normalsSM = getNormals1Face(1:nPoints,template.faces,normalsSM);

%% compute NN
[idxsNN, distNN] = knnsearch(scan.points,pointsSM);

%% find scan normals
assert(size(idxsNN,1) == nPoints);
normalsScanAll = getNormals1Face(scan.pointsIdxs,scan.faces,scan.normals);
normalsScan = normalsScanAll(idxsNN,:);

%% check the angle between normals
isValidNN = checkAngle(normalsScan,normalsSM,threshNormAngle);
    
% do not register open to closed hands
isValidNN(idxHand) = 0;

if (~isempty(template.idxsUse))
    % use only subset of vertices
    isValidNN = isValidNN.*template.idxsUse;
end

distAll = distNN' * isValidNN;
err = distAll/sum(isValidNN);

%% optimization parameters
[options,poseLB,poseUB,shapeLB,shapeUB] = getOptionsOptimizer(modelDir,template);

poseLB(template.poseParamsIgnoreIdxs) = 0;
poseUB(template.poseParamsIgnoreIdxs) = 0;

if (~bScale)
    poseLB(end) = template.poseParams(end);
    poseUB(end) = template.poseParams(end);
    shapeLB(end) = template.poseParams(end);
    shapeUB(end) = template.poseParams(end);
end

eps_err = 1e-3;
errPrev = err+eps_err+1;

%% optimization loop
while abs(errPrev - err) > eps_err    
    
    errPrev = err;
    
    fprintf('fit pose\n');
    [poseParams, ~] = fmincon(@PoseFunc,poseParams,[],[],[],[],poseLB,poseUB,[],options);
    [pointsSM, ~] = shapepose(poseParams(1:end-1),shapeParams(1:end-1),evectors,modelDir);
    sc = poseParams(end);
    pointsSM = sc * pointsSM;
    idxsNN = knnsearch(scan.points,pointsSM);
    
    % check the angle between normals
    normalsScan = normalsScanAll(idxsNN,:);
    normalsSM = getNormals(pointsSM, template.faces);
    normalsSM = getNormals1Face(1:nPoints,template.faces,normalsSM);
    isValidNN = checkAngle(normalsScan,normalsSM,threshNormAngle);
    
    % do not register open to closed hands
    isValidNN(idxHand) = 0; 
   
    if (~isempty(template.idxsUse))
        % use only subset of vertices
        isValidNN = isValidNN.*template.idxsUse;
    end
    fprintf('sum(isValidNN): %1.1f\n',sum(isValidNN));
    if (bFitShape)
        fprintf('fit shape\n');
        shapeParams(end) = sc;
        [shapeParams, ~] = fmincon(@ShapeFunc, shapeParams,[],[],[],[],shapeLB,shapeUB,[],options);
        % new model code
        [pointsSM, ~] = shapepose(poseParams(1:end-1),shapeParams(1:end-1),evectors,modelDir);
        sc = shapeParams(end);
        pointsSM = sc * pointsSM;
        [idxsNN, distNN] = knnsearch(scan.points,pointsSM);
        
        % check the angle between normals
        normalsScan = normalsScanAll(idxsNN,:);
        normalsSM = getNormals(pointsSM,template.faces);
        normalsSM = getNormals1Face(1:nPoints,template.faces,normalsSM);
        isValidNN = checkAngle(normalsScan,normalsSM,threshNormAngle);
        
        % do not register open to closed hands
        isValidNN(idxHand) = 0;
        
        if (~isempty(template.idxsUse))
            % use only subset of vertices
            isValidNN = isValidNN.*template.idxsUse;
        end
    end
    
    distAll = distNN' * isValidNN;
    err = distAll / sum(isValidNN);
    
    poseParams(end) = sc;
end

%% output
pointsSM_best = shapepose(poseParams(1:end-1), shapeParams(1:end-1),evectors,modelDir);

sc = shapeParams(end);
pointsSM_best = sc * pointsSM_best;
template.shapeParams = shapeParams(1:end-1);
template.poseParams = poseParams;
template.points = pointsSM_best;
template.pointsIdxsScanNN = idxsNN;
template.pointsHasValidNN = isValidNN;
fprintf('done\n');

    %% fit pose objective function
    function E = PoseFunc(poseParam)
        scTmp = poseParam(1,32);
        [pointsSM, ~] = shapepose(poseParam(1:end-1), shapeParams(1:end-1),evectors,modelDir);
        pointsSM = scTmp * pointsSM;
        pointsScanNN = scan.points(idxsNN, :);
        dist = sqrt(sum((pointsScanNN - pointsSM) .^ 2, 2));
        E = sum(dist.*isValidNN);
    end

    %% fit shape objective function
    function E = ShapeFunc(shapeParam)
        scTmp = shapeParam(end);
        [pointsSM, ~] = shapepose(poseParams(1:end-1), shapeParam(1:end-1),evectors,modelDir);
        pointsSM = scTmp * pointsSM;
        pointsScanNN = scan.points(idxsNN, :);
        dist = sqrt(sum((pointsScanNN - pointsSM) .^ 2, 2));
        E = sum(dist.*isValidNN);
    end
end
