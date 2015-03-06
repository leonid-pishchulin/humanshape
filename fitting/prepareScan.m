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

function scan = prepareScan(scan,scanFilenames)

rng('default');
rng(42);

for i = 1:length(scanFilenames)
    %% read scan
    [~,~,ext] = fileparts(scanFilenames{1});
    if (strcmp(ext,'.obj'))
        [pointsAllScan,facesScan] = read_obj(scanFilenames{i});
        confidenceScan = ones(size(pointsAllScan,1),1);
    elseif (strcmp(ext,'.mat'))
        load(scanFilenames{i});
        assert(exist('points','var')>0);
        pointsAllScan = points;
        if (exist('faces','var') == 0)
            faces = zeros(0,3);
        end
        facesScan = faces;
        confidenceScan = ones(size(pointsAllScan,1),1);
    else
        [pointsAllScan,confidenceScan,facesScan] = read_ply(scanFilenames{i});
    end
    
    pointsAllScan = m2mm(pointsAllScan);
    
    %% subsample scan
%     nPointsSample = 16000;
    nPointsSample = min(6449*3,size(pointsAllScan,1));
    fprintf('nPointsSample: %d/%d\n',nPointsSample,size(pointsAllScan,1));
    
%     if (i == 1)
        pointsIdxsScan = randperm(size(pointsAllScan,1),nPointsSample);
%     end
    
    points = pointsAllScan(pointsIdxsScan(1:nPointsSample),:);
    
    %% compute normals
    normalsScan = getNormals(pointsAllScan,facesScan);
    
    scan(i).points = points;
    scan(i).normals = normalsScan;
    scan(i).faces = facesScan;
    scan(i).confidence = confidenceScan;
    scan(i).pointsIdxs = pointsIdxsScan;
end

end