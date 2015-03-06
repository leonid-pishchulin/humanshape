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

function [options,poseLB,poseUB,shapeLB,shapeUB] = getOptionsOptimizer(modelDir,template)

options = optimset('Display', 'iter');
options.MaxFunEvals = 10000;
options.DiffMinChange = 0.0001;
options.Algorithm = 'interior-point';
options.LargeScale = 'on';

if (isfield(template, 'poseLB') && isfield(template, 'poseUB') && ~isempty(template.poseLB) && ~isempty(template.poseUB))
    poseUB = template.poseUB;
    poseLB = template.poseLB;
else
    % global and local rotations
    poseLB = ones(1,length(template.poseParams))*(-2*pi);
    poseUB = ones(1,length(template.poseParams))*2*pi;
    
    % translation
    poseLB(4:6) = -inf;
    poseUB(4:6) = inf;
    
    % scale
    poseLB(length(template.poseParams)) = 0;
    poseUB(length(template.poseParams)) = inf;
end
% space of human shapes: meanShape +/- 3 sigma
if (template.nPCA <= 20)
    load([modelDir '/evalues'], 'evalues');
else
    load([modelDir '/evalues-all'], 'evalues');
end
assert(template.nPCA <= length(evalues));

shapeLB = -3*sqrt(evalues(1:length(template.shapeParams)));
shapeUB =  3*sqrt(evalues(1:length(template.shapeParams)));

% scale
shapeLB(length(template.shapeParams)+1) = 0;
shapeUB(length(template.shapeParams)+1) = inf;

end