function [idxAbove, idxBelow] = getHeadTopIdx(templatePoints,points)

%% hardcode the cutting plane by using 3 points
idx1 = 4122;
idx2 = 4082;
idx3 = 982;

p1 = templatePoints(idx1,:);
p2 = templatePoints(idx2,:);
p3 = templatePoints(idx3,:);

% Getting the Normal to the cutting plane & the constant term of the cutting plane equation
normal = cross(p1-p2, p1-p3);
normal = normal./norm(normal);
offset = -1* (p1 * transpose(normal));

nPoints = size(points,1);
% split SM points into above cutting plane and below cutting plane
idxAbove = [];
for i=1:nPoints
    if ((points(i,:) * transpose(normal)) + offset < 0)
%         if ((points(i,:) * transpose(normal)) + offset >= 0)
        idxAbove = cat(2,idxAbove,i);
    end
end
i=1:1:nPoints;
idxBelow = setdiff(i,idxAbove);
    
end