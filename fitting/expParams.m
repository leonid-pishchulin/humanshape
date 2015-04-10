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

function p = expParams(expidx)
p = [];

% root directory of code
p.rootDir = '/home/leonid/code/humanshape/';

% landmarks used to fit the model
p.landmarksSM = {[p.rootDir 'fitting/landmarksIdxs73.mat']};
% faces for model visualization
p.facesSM = [p.rootDir 'fitting/facesShapeModel.mat'];
% model mesh with skining weights
p.modelMesh = [p.rootDir 'shapemodel/model.dat'];
% root directory for experiments
p.expDir = [p.rootDir 'experiments/'];
% directory to save visualizations
p.figures = '';

switch expidx
    
    case 0
        p.name = 'init-nrd-3';
        % directory with fitting results
        % e.g. p.fitDir = 'caesar-norm-nh-fitted-meshes';
        % each fitted mesh should be stored as p.fitDir/<mesh-name>/NRD.mat
        p.fitDir = [p.expDir '/' p.name];
        %load initial fitting, if available
        p.initDir = p.fitDir;
        % directory to load existing PCA model from
        p.modelInDir = [p.expDir 'models/caesar'];
        % directory to save a newly learned PCA model to
        p.modelOutDir = [p.fitDir '/model'];
        % shape and pose fitting prior to non-rigid deformation (NRD)
        p.bInit = true;
        % variants of NRD:
        % -1 - no NRD; 0 - NRD, const weights; 1-3 - NRD, reduce weights
        % p.nrdWidx = 3 leads to best results
        p.nrdWidx = 3;
        % number of shape space parameters
        p.nPCA = 20;
        % color id for visualization
        p.colorIdxs = [6 1];
end

if (isfield(p,'colorIdxs') && ~isempty(p.colorIdxs))
    p.colorName = eval_get_color(p.colorIdxs);
    p.colorName = p.colorName ./ 255;
end

end