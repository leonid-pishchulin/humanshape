function visLandmarks(scan,template)

% load([modelDir '/evectors.mat'],'evectors');
% template.points = changeShapePose(template.poseParams(end-1),template.shapeParams,evectors,modelDir);

labels = cell(size(scan.landmarks,1),1);
for i = 1:length(labels)
    labels{i} = num2str(i);
end

landmarksScan = scan.landmarks(scan.landmarkIdxs,:);
landmarksTemplate = template.points(template.landmarksIdxs(scan.landmarkIdxs),:);
labels = labels(scan.landmarkIdxs);

clf;
plot3(landmarksScan(:,1),landmarksScan(:,2),landmarksScan(:,3),'g+','markerSize',12); hold on;
text(landmarksScan(:,1),landmarksScan(:,2),landmarksScan(:,3),labels); hold on;
plot3(landmarksTemplate(:,1),landmarksTemplate(:,2),landmarksTemplate(:,3),'r*','markerSize',12); hold on;
text(landmarksTemplate(:,1),landmarksTemplate(:,2),landmarksTemplate(:,3),labels);
axis equal;
grid on;
view(0,90);

end