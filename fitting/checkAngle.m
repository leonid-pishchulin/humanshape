function res = checkAngle(vec1,vec2,threshAngle)

if (isempty(vec1) || isempty(vec2))
    res = ones(size(vec2,1),1);
    return;
end

crossProd = cross(vec1,vec2,2);
dotProd = sum((vec1.*vec2),2);
angle = atan2((sum((crossProd.^2),2)).^0.5,dotProd)/pi*180;

for i = 1:length(angle)
    while(angle(i) < 0)
        angle(i) = angle(i) + 360;
    end
    
    while(angle(i) > 360)
        angle(i) = angle(i) - 360;
    end
end
res = (angle <= threshAngle);
% fprintf('checkAngle(); thresh: %f, checkPassed: %d/%d\n', thresh, sum(checkPassed), length(angle));
end