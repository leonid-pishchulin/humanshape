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

function fitMesh(scanFilenames,landmarksScanFilename,expidx)
fprintf('fitMesh()\n');

% scanFilenames = {'scan.ply'};
% landmarksFilenames = {'scan.lnd'};
% expidx = 1;

if (ischar(expidx))
    expidx = str2double(expidx);
end

p = expParams(expidx);

if (~exist(p.fitDir,'dir'))
    mkdir(p.fitDir);
end

% remove mosek from matlab path, if any
s = which('quadprog.m');
if ~isempty(strfind(s,'mosek')) && ~isdeployed
    new_path = removeMosek();
    matlabpath(new_path);
end

% create directories for fitting results
saveDir = cell(length(scanFilenames),1);
initDir = cell(length(scanFilenames),1);
for i = 1:length(scanFilenames)
    [~,nameScan] = fileparts(scanFilenames{i});
    saveDir{i} = [p.fitDir '/' nameScan];
    if (~exist(saveDir{i}, 'dir'))
        mkdir(saveDir{i});
    end
    initDir{i} = [p.initDir '/' nameScan];
    assert(exist(initDir{i},'dir')>0);
end

% print parameters
fprintf('scanFilename: %s\n',scanFilenames{1});
fprintf('landmarksFilename: %s\n',landmarksScanFilename{1});
fprintf('saveDir: %s\n',saveDir{1});
fprintf('modelDir: %s\n',p.modelInDir);
fprintf('bInit: %d\n',p.bInit);
fprintf('initDir: %s\n',p.initDir);
fprintf('nPCA: %d\n',p.nPCA);

% read landmarks
[scan,template] = readLandMarksAll(landmarksScanFilename,p.landmarksSM);

% prepare template
load(p.facesSM,'faces');
for i = 1:length(template)
    template(i).faces = faces;
    
    DOF = 30;
    template(i).poseParams = zeros(1,DOF+2);
    template(i).poseParams(end) = 1;
    load([p.modelInDir '/evectors'], 'evectors');
    
    assert(p.nPCA <= size(evectors,1));
    template(i).shapeParams = zeros(1,p.nPCA);
    load([p.modelInDir '/meanShape'], 'points');
    % center meanShape
    points = points - repmat(mean(points),size(points,1),1);
    template(i).points = points;
    template(i).nPCA = p.nPCA;
end

% prepare scan
scan = prepareScan(scan,scanFilenames);
assert(length(scan) == length(template));

for i=1:length(scan)
    % save subsampled scan
    points = scan(i).points;
    faces = scan(i).faces;
    template(i).idxsUse = [];
    template(i).poseParamsIgnoreIdxs = [23 28]; % ignore wrist rotations
    normals = getNormals1Face(scan(i).pointsIdxs,scan(i).faces,scan(i).normals);
    save([saveDir{i} '/scan'],'points','normals');
end

if (length(scan) == 1)
    % register the scan
    template = register(scan,template,saveDir{1},p.modelInDir,p.bInit,p.nrdWidx,initDir{1});
else
    % joint optimization over shape parameters when multiple scans available
    template = registerJoint(scan,template,saveDir,modelDir,initDir);
end    
close all;

end
