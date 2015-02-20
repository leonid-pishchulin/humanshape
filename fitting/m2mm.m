function points = m2mm(points)

maxDist = max(max(points) - min(points));
% make sure the height is mm
% TODO: better way?
if (maxDist < 500)
    points = points * 1000;
end

end