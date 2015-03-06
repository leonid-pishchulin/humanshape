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