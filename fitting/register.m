function template = register(scan,template,saveDir,modelDir,bInit,nrdWidx,initDir)
fprintf('\nregister()\n');

if (nargin < 7)
    initDir = saveDir;
end

% rigid alignment of template to the scan
scan = rigidAlignTemplate2Scan(scan,template);
% visualize landmarks and fit
visLandmarks(scan,template);
visFit(scan,template);

if (bInit)
    % fit model shape and pose parameters
    template = fitModel(scan,template,modelDir,initDir);
end

if (nrdWidx > -1)
    % perform non-rigid deformation (NRD)
    fname = [saveDir '/NRD'];
    try
        fprintf('load %s\n', fname);
        load(fname,'points');
    catch
        tic
        [template,dist] = NRD(scan,template,nrdWidx);
        points4     = [template.points ones(size(template.points,1),1)]/scan.T';
        points      = points4(:,1:3);
        Amatrix     = template.Amatrix;
        pointsIdxsScanNN = template.pointsIdxsScanNN;
        pointsHasValidNN = template.pointsHasValidNN;
        T = scan.T;
        save(fname,'points','Amatrix','pointsIdxsScanNN','pointsHasValidNN','dist','T');
        toc
    end
end
end