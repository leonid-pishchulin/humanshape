%{
    This file is part of the evaluation of the 3D human shape model as described in the paper:

    Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobalt and Bernt Schiele
    Building Statistical Shape Spaces for 3D Human Modeling
    ArXiv, March 2015

    Please cite the paper if you are using this code in your work.
    
    Contributors: Stefanie Wuhrer, Leonid Pishchulin.

    The code may be used free of charge for non-commercial and
    educational purposes, the only requirement is that this text is
    preserved within the derivative work. For any other purpose you
    must contact the authors for permission. This code may not be
    redistributed without permission from the authors.
%}

function learnPCA(expidx)

fprintf('learnPCA()\n');

p = expParams(expidx);
assert(exist(p.fitDir, 'dir') > 0);

if (~exist(p.modelOutDir, 'dir'))
    mkdir(p.modelOutDir);
end

% read fittings
fnameData = [p.fitDir '/alignedScans.mat'];
if (exist(fnameData, 'file') > 0)
    load(fnameData,'data');
else
    data = readFits(p.fitDir);
    if (size(data,1) < 2)
        warning(['At least two fits required, found: ' num2str(size(data,1)) '. Exiting']);
        return;
    end
    tic
    % procrustes analysis to rigidly align the scans
    data = rigidAlign(data',size(data,2),size(data,1))';
    toc
    fprintf('save: %s\n',fnameData);
    save(fnameData,'data');
end

totalSamples = size(data,1);
dimensions = size(data,2);

% learn PCA model
[~, evalues, evectors, meanData] = ErrorEvaluation(0,data',[],dimensions,totalSamples,0,totalSamples-1,totalSamples,0);
if (~exist(p.modelOutDir, 'dir'))
    mkdir(p.modelOutDir);
end
save([p.modelOutDir '/evectors-flat.mat'], 'evectors');
save([p.modelOutDir '/meanData-flat.mat'], 'meanData');
save([p.modelOutDir '/evalues.mat'], 'evalues');

fprintf('done\n');
end
