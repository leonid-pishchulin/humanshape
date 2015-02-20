function [template, dist] = fitPose(scan,template,modelDir,bScale)
fprintf('\nfitPose()\n');

if (nargin < 4)
    bScale = true;
end

landmarksScan = scan.landmarks(scan.landmarkIdxs,:);
landmarksIdxsSM = template.landmarksIdxs(scan.landmarkIdxs);

%% optimization options
[options,poseLB,poseUB] = getOptionsOptimizer(modelDir,template);

poseLB(template.poseParamsIgnoreIdxs) = 0;
poseUB(template.poseParamsIgnoreIdxs) = 0;

if (~bScale)
    poseLB(end) = template.poseParams(end);
    poseUB(end) = template.poseParams(end);
end

x0 = template.poseParams; % model parameteres
poseParams = template.poseParams;

% mean shape
shapeParams = template.shapeParams;

% load eigenvectors
load([modelDir '/evectors'], 'evectors');
assert(template.nPCA <= size(evectors,1));
evectors = evectors(1:template.nPCA,:);
% visLandmarks(scan,template);

%% optimization
optPoseParam = fmincon(@optFunc,x0,[],[],[],[],poseLB,poseUB,[],options);
x0 = optPoseParam;
pointsSM = shapepose(optPoseParam(1:end-1),shapeParams,evectors,modelDir);
pointsSM = pointsSM * x0(1,end);

%% output
[pointsSM_best, ~] = shapepose(optPoseParam,shapeParams,evectors,modelDir);
pointsSM_best = pointsSM_best * optPoseParam(1,32);
template.points = pointsSM_best;
template.poseParams = optPoseParam;
dist = sqrt(sum((landmarksScan - pointsSM_best(landmarksIdxsSM,:)).^2,2));

fprintf('done\n');

    %% objective function
    function E = optFunc(params)
        poseParams = params(1:end-1);
        [pointsSM, ~] = shapepose(poseParams,shapeParams,evectors,modelDir);
        sc = params(1,end);
        pointsSM = sc * pointsSM;
        E = sum(sqrt(sum((landmarksScan - pointsSM(landmarksIdxsSM,:)).^2,2)));
    end
end