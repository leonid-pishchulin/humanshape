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