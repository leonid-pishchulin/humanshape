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