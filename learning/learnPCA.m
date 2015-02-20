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
