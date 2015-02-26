expidx = 0;
p = expParams(expidx);

%% manipulate shape and pose
fprintf('mean pose and shape, press any key\n');
poseParams = zeros(1,31); % mean pose
shapeParams = zeros(1,p.nPCA); % mean shape
load([p.modelInDir '/evectors'],'evectors'); % shape space eigenvectors
load([p.modelInDir '/evalues'],'evalues'); % shape space eigenvalues
evectors = evectors(1:p.nPCA,:);
evalues = evalues(1:p.nPCA);
points = shapepose(poseParams,shapeParams,evectors,p.modelInDir);

% show model
load(p.facesSM,'faces');
clf;
showmodel(points,faces,'r',[],0);
axis equal; view(45,22.5); pause; 

fprintf('change pose and shape, press any key\n');
% description of pose parameters in shapemodel/poseParamsDescript.m
poseParams(24) = 45/180*pi; % rotate right shoulder by 45 degrees
shapeParams(1) = 3*sqrt(evalues(1)); % change height to 3 st.d. from mean
points = shapepose(poseParams,shapeParams,evectors,p.modelInDir);

% show model
hold on;
showmodel(points,faces,'g',[],0);
axis equal; view(45,22.5); pause;

fprintf('visualize eigenvector scaled by 3 st.d., press any key\n');
vpAngle = 90; % visualization viewpoint
idxShape = 1; % first eigenvector
sign = -1; % sign of st.d.
bSave = false; % save figure
clf;
visModelVP(expidx,vpAngle,idxShape,sign,bSave); pause;

%% register human scan
scanName = 'scan';
scanFilenames = {[scanName '.ply']};
landmarkFilenames = {[scanName '.lnd']};
fprintf('register human scan\n');
fitMesh(scanFilenames,landmarkFilenames,expidx);

% visualize registration
fprintf('results: red - fitted mesh, blue - overlaid scan\n');
fprintf('result of pose fitting using landmarks, press any key\n');
clf; visFitDir([p.fitDir '/' scanName],0); pause;
fprintf('result of pose and shape fitting using all vertices, press any key\n');
clf; visFitDir([p.fitDir '/' scanName],1); pause;
fprintf('result of non-rigid deformation, press any key\n');
clf; visFitDir([p.fitDir '/' scanName],2); pause;

%% learn PCA model
% learnPCA(expidx);
% createModel(expidx);
fprintf('done\n');