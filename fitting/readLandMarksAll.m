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

function [scan,template,bIsBad] = readLandMarksAll(landmarksScanFilenames,landmarksSMFilenames)

[~,~,ext] = fileparts(landmarksScanFilenames{1});
if (strcmp(ext,'.mat'))
    load(landmarksScanFilenames{1},'landmarks');
else
    landmarks = readLandmarks(landmarksScanFilenames{1});
end

nLandmarks = min(size(landmarks,1),73);
landmarks = landmarks(1:nLandmarks,:);

assert(length(landmarksScanFilenames) == length(landmarksSMFilenames));

% CAESAR: use first 72 landmarks
if nLandmarks == 73
    nLandmarksUse = 72;
else
    nLandmarksUse = nLandmarks;
end

bIsBad = true(nLandmarks,1);
bIsBad(1:nLandmarksUse) = false;

scan.landmarks = [];

for i = 1:length(landmarksScanFilenames)
    if (strcmp(ext,'.mat'))
        load(landmarksScanFilenames{i},'landmarks');
    else
        landmarks = readLandmarks(landmarksScanFilenames{i});
    end
    
    load(landmarksSMFilenames{i},'landmarksIdxs');
    template(i).landmarksIdxs = landmarksIdxs;

    landmarks = landmarks(1:nLandmarks,:);
    landmarks = m2mm(landmarks);
    scan(i).landmarks = landmarks;
    scan(i).landmarkIdxs = find(~isnan(scan(i).landmarks(:,1)) & ~isnan(template(i).landmarksIdxs));
end

end