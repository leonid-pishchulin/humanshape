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

function [template, dist] = NRD(scan,template,nrdWidx)

fprintf('\nNRD()\n');

if (nargin < 3)
    nrdWidx = 3;
end

global NUM_LM;
global NUM_TL_POINTS;
global NUM_SC_POINTS;

global LM_TL_IDX;
global NN_VERT_TL;

global TL_POINTS;
global LM_SC;
global LM_TL;

global NN_SC;
global CONF_SC_POINTS;

global TL_HEAD_IDX;
global TL_HAND_IDX;
global TL_ALL_IDX;

global NUM_TL_HEAD;
global NUM_TL_HAND;
global NUM_TL_ALL;

global IS_VALID_NN;

global W_DATA;
global W_SMOOTH_HEAD;
global W_SMOOTH_HANDS;
global W_SMOOTH_GEN;
global W_LM;

global E_TOTAL;
global E_DATA;
global E_SMOOTH;
global E_LM;

%% separate landmarks for the head and hands
LM_TL_IDX = template.landmarksIdxs(scan.landmarkIdxs);
NUM_SC_POINTS = size(scan.points,1);

NUM_LM = size(scan.landmarks(scan.landmarkIdxs,:),1);
NUM_TL_POINTS = size(template.points,1);

% head and hand vertex idxs
load('VertexIdxSpecParts','idxHead','idxHand','idxFoot');
TL_HEAD_IDX = idxHead;
TL_HAND_IDX = idxHand;

TL_ALL_IDX = 1:1:NUM_TL_POINTS;
TL_ALL_IDX = setdiff(TL_ALL_IDX,TL_HEAD_IDX);
TL_ALL_IDX = setdiff(TL_ALL_IDX,TL_HAND_IDX);

NUM_TL_HEAD = size(TL_HEAD_IDX,2);
NUM_TL_HAND = size(TL_HAND_IDX,2);
NUM_TL_ALL = size(TL_ALL_IDX,2);

%% init variables for NRD
TL_POINTS = [];
LM_SC = [];
LM_TL = [];
Amatrix = []; % Deformation matrix
lowerbound = []; % lowerbound for Amatrix optimization
upperbound = []; % upperbound for Amatrix optimization

for i=1:NUM_TL_POINTS
    TL_POINTS = cat(3,TL_POINTS,[template.points(i,:) 1]');
end

for i = 1:NUM_LM
    LM_SC = cat(3,LM_SC,[scan.landmarks(scan.landmarkIdxs(i),:) 1]');
    LM_TL = cat(3,LM_TL,[template.points(LM_TL_IDX(i,1),:) 1]');
end

for i = 1:NUM_TL_POINTS
    Amatrix = cat(1,Amatrix,reshape(eye(4,4),[16 1]));
    lowerbound = cat(1,lowerbound,-10000 * ones(16,1));
    upperbound = cat(1,upperbound,10000 * ones(16,1));
end

scan.normals1face = getNormals1Face(scan.pointsIdxs,scan.faces,scan.normals);
scan.pointsIdxsFaces = find(~isnan(scan.normals1face(:,1)));

normalsTemplate = getNormals(template.points,template.faces);
normalsTemplate = getNormals1Face(1:NUM_TL_POINTS,template.faces,normalsTemplate);
template.normals = normalsTemplate;

NN_VERT_TL = getVertexNN(template.faces);

pointsIdxsScanNN = getNNheadAndBody(scan,template);

IS_VALID_NN = ones(NUM_TL_POINTS,1);
IS_VALID_NN(isnan(pointsIdxsScanNN)) = 0;
IS_VALID_NN(TL_HAND_IDX) = 0;
IS_VALID_NN(idxFoot) = 0;

pointsIdxsScanNN(isnan(pointsIdxsScanNN)) = 1;

NN_SC = scan.points(pointsIdxsScanNN,:);
CONF_SC_POINTS = scan.confidence(pointsIdxsScanNN,:);

NN_SC = NN_SC .* repmat(IS_VALID_NN, 1, 3);
CONF_SC_POINTS = CONF_SC_POINTS .* IS_VALID_NN;

[wrSmoothBody,wrSmoothHand,wrSmoothHead,wrLandmarks] = getWeightReduceFactor(nrdWidx);
fprintf('nrdWidx: %d\n',nrdWidx);

%% initial NRD: LM + SMOOTH are considered 
W_DATA = 0;
W_SMOOTH_HEAD = 1e+6;
W_SMOOTH_HANDS = 1e+6;
W_SMOOTH_GEN = 1e+6;
W_LM = 1e-3;

pointsTemplate = template.points;

fprintf('W_DATA: %1.2f; W_SMOOTH_HEAD: %1.2f; W_SMOOTH_HANDS: %1.2f; W_SMOOTH_GEN: %1.2f; W_LM: %1.2f\n', W_DATA, W_SMOOTH_HEAD, W_SMOOTH_HANDS, W_SMOOTH_GEN, W_LM);
distInit = sqrt(sum((pointsTemplate - NN_SC).^2,2))'*IS_VALID_NN;
fprintf('distInit: %f \n', distInit);
fprintf('init; ');

%which lbfgsb

try
    Amatrix = lbfgsb(Amatrix,lowerbound,upperbound,'NRDFunction','NRDgradientFunction');
catch exception
    disp(exception.message);
end

Amatrix = reshape(Amatrix,[4 4 NUM_TL_POINTS]);

dist = zeros(NUM_TL_POINTS,1);
for i=1:NUM_TL_POINTS
    smp = Amatrix(:,:,i)*TL_POINTS(:,:,i);
    pointsTemplate(i,:) = smp(1:3);
    dist(i) = norm(Amatrix(:,:,i)*TL_POINTS(:,:,i) - [NN_SC(i,:)';1]);
end
distAll = dist'*IS_VALID_NN;

normalsTemplate = getNormals(pointsTemplate,template.faces);
normalsTemplate = getNormals1Face(1:NUM_TL_POINTS,template.faces,normalsTemplate);
template.normals = normalsTemplate;

err = distAll/sum(IS_VALID_NN);

Amatrix = reshape(Amatrix,[16*NUM_TL_POINTS 1]);
fprintf('ERR: %1.2f; NUM_VALID_NN: %d; distAll: %1.2f; E_DATA: %1.2f; E_SMOOTH: %1.2f; E_LM: %1.2f; E_TOTAL: %1.2f \n', distAll/sum(IS_VALID_NN), sum(IS_VALID_NN), distAll, E_DATA, E_SMOOTH, E_LM, E_TOTAL);

%% NRD: DATA + LM + SMOOTH 
eps_err = 1e-3;
errPrev = err+eps_err+1;
W_DATA = 1;
W_SMOOTH_HEAD = 1e+6;
W_SMOOTH_HANDS = 1e+6;
W_SMOOTH_GEN = 1e+6;
W_LM = 1e-3;

iter = 0;
fprintf('W_DATA: %1.2f; W_SMOOTH_HEAD: %1.2f; W_SMOOTH_HANDS: %1.2f; W_SMOOTH_GEN: %1.2f; W_LM: %1.2f\n', W_DATA, W_SMOOTH_HEAD, W_SMOOTH_HANDS, W_SMOOTH_GEN, W_LM);

while (abs(errPrev - err) > eps_err && W_SMOOTH_GEN >= W_DATA*1e+3)
    
    iter = iter + 1;
    fprintf('iter: %d; ', iter);
    
    errPrev = err;
    
    % NN
    template.points = pointsTemplate;
    template.normals = normalsTemplate;
    pointsIdxsScanNN = getNNheadAndBody(scan,template);
    IS_VALID_NN = ones(NUM_TL_POINTS,1);
    IS_VALID_NN(isnan(pointsIdxsScanNN)) = 0;
    IS_VALID_NN(TL_HAND_IDX) = 0;
    IS_VALID_NN(idxFoot) = 0;
    
    pointsIdxsScanNN(isnan(pointsIdxsScanNN)) = 1;
    
    NN_SC = scan.points(pointsIdxsScanNN,:); % nearest neighbours of template points on scan
    CONF_SC_POINTS = scan.confidence(pointsIdxsScanNN,:); 

    try
        Amatrix = lbfgsb(Amatrix,lowerbound,upperbound,'NRDFunction','NRDgradientFunction');
    catch exception
        disp(exception.message);
    end
    
    Amatrix = reshape(Amatrix,[4 4 NUM_TL_POINTS]);
    
    dist = zeros(NUM_TL_POINTS,1);
    for i=1:NUM_TL_POINTS
        smp = Amatrix(:,:,i)*TL_POINTS(:,:,i);
        pointsTemplate(i,:) = smp(1:3);
        dist(i) = norm(Amatrix(:,:,i)*TL_POINTS(:,:,i) - [NN_SC(i,:)';1]);
    end
    
    distAll = dist'*IS_VALID_NN;
    
    normalsTemplate = getNormals(pointsTemplate,template.faces);
    normalsTemplate = getNormals1Face(1:NUM_TL_POINTS,template.faces,normalsTemplate);
    
    err = distAll/sum(IS_VALID_NN);
    Amatrix = reshape(Amatrix,[16*NUM_TL_POINTS 1]);
    
    % reduce weights
    W_SMOOTH_GEN   = W_SMOOTH_GEN   * wrSmoothBody;
    W_SMOOTH_HANDS = W_SMOOTH_HANDS * wrSmoothHand;
    W_SMOOTH_HEAD  = W_SMOOTH_HEAD  * wrSmoothHead;
    W_LM = W_LM * wrLandmarks;
    
    fprintf('ERR: %1.2f; NUM_VALID_NN: %d; AggDist: %1.2f; E_DATA: %1.2f; E_SMOOTH: %1.2f; E_LM: %1.2f; E_TOTAL: %1.2f \n', err, sum(IS_VALID_NN), distAll, E_DATA, E_SMOOTH, E_LM, E_TOTAL);
    
end
Amatrix = reshape(Amatrix,[4 4 NUM_TL_POINTS]);
template.Amatrix = Amatrix;
template.points = pointsTemplate;
[pointsIdxsScanNN, dist] = knnsearch(scan.points,template.points);
template.pointsIdxsScanNN = pointsIdxsScanNN;
template.pointsHasValidNN = IS_VALID_NN;
end