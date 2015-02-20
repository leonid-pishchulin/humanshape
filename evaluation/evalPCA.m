function evalPCA(evalMode,id,saveDir,registerDir,nSamples,nTraininMeshes,bReadNils,dbDir,scans_info_filename,bInitOnly)

fprintf('evalPCA()\n');
% input Data
% evalMode    1-compactness 2-generalization 3-specificity  4-learnDistribution
% id            generalization: leave-one-out id; specificity: num of PCA
% saveToDir     Path to save mat files
% dbDir       = '/scratch/BS/pool0/monica/CAESAR-DB/CAESAR/'
% registerDir = '/BS/leonid-people-3d/work/experiments-caesar/register3D-iter2/'
% goodResDir  = '/BS/leonid-people-3d/work/experiments-caesar/register3D-iter2/vis_all'

if (~isdeployed)
    addpath('./errorEval/');
end

if ischar(evalMode)
    evalMode = str2num(evalMode);
end

if ischar(id)
    id = str2num(id);
end

if nargin < 5
    nSamples = inf;
elseif (ischar(nSamples))
    nSamples = str2num(nSamples);
end

if nargin < 6
    nTraininMeshes = inf;
elseif (ischar(nTraininMeshes))
    nTraininMeshes = str2num(nTraininMeshes);
end

if nargin < 7
    bReadNils = false;
end

if nargin < 8
    dbDir = '/BS/databases00/CAESAR-CT';
end

if nargin < 9
%     scans_info_filename = '/BS/leonid-people-3d/work/experiments-caesar/data/scans_info.mat';
    scans_info_filename = '/BS/leonid-people-3d/work/experiments-caesar/data/scans.mat';
end

if nargin < 10
    bInitOnly = false;
end

fprintf('evalMode: %d\n', evalMode);
fprintf('id: %d\n', id);
fprintf('saveToDir: %s\n', saveDir);
fprintf('registerDir: %s\n', registerDir);
fprintf('nSamples: %d\n', nSamples);
fprintf('nTraininMeshes: %d\n', nTraininMeshes);
fprintf('dbDir: %s\n', dbDir);
fprintf('scans_info_filename: %s\n', scans_info_filename);
fprintf('bInitOnly: %d\n',bInitOnly);

assert(exist(dbDir, 'dir') > 0);
assert(exist(registerDir, 'dir') > 0);

% not used in the code 
numComponents = 0;

if (~exist(saveDir, 'dir'))
    mkdir(saveDir);
end

load(scans_info_filename, 'SCANS');

% read scans
fnameData = [saveDir '/nrdAll.mat'];
if (exist(fnameData, 'file') > 0)
    load(fnameData,'data','scan_idxs');
else
    if (~bReadNils)
        [data,scan_idxs] = readNRD(registerDir,SCANS,bInitOnly);
        if (nSamples < length(scan_idxs))
            [male_idxs,female_idxs] = sampleData(nSamples,SCANS,scan_idxs);
            idxs = sort([male_idxs;female_idxs]);
            idxsRel = zeros(length(idxs),1);
            for i=1:length(idxsRel)
                idxsRel(i) = find(scan_idxs == idxs(i));
            end
            data = data(idxsRel,:);
            %         if (nSamples <= 100)
            % visualize samples
            %             visualizeScans(data,[saveDir '/vis'],SCANS.name(idxs));
            %         end
        end
    else
        bSampleSimple = true;
        [data,scan_idxs] = readNRDnils(registerDir);
        if (nSamples < length(scan_idxs) && bSampleSimple)
            [idxs] = sampleDataSimple(nSamples,scan_idxs);
            data = data(idxs,:);
        elseif (nSamples < length(scan_idxs) && ~bSampleSimple)
            [male_idxs,female_idxs] = sampleData(nSamples,SCANS,scan_idxs);
            idxs = sort([male_idxs;female_idxs]);
            idxsRel = zeros(length(idxs),1);
            for i=1:length(idxsRel)
                idxsRel(i) = find(scan_idxs == idxs(i));
            end
            data = data(idxsRel,:);
        end
        
    end
    
    tic
    data = rigidAlign(data',size(data,2),size(data,1))';
    toc
    fprintf('save: %s\n',fnameData);
    save(fnameData,'data','scan_idxs');
end

if (evalMode == 4 || evalMode == 5)
    fnameScans = '/BS/leonid-pose/work/experiments-caesar/model-mpii-all-samples-warm-sched-3-bootstrap-4/scansAll.mat'; % always use scans after bootstrap 4
    if (exist(fnameScans, 'file') > 0)
        load(fnameScans,'dataGT','scan_idxs_gt');
        assert(sum(scan_idxs_gt-scan_idxs) == 0);
    else
        [dataGT,scan_idxs_gt] = readScans(registerDir,SCANS);
        assert(sum(scan_idxs_gt-scan_idxs) == 0);
        tic
        dataGT = rigidAlign(dataGT',size(dataGT,2),size(dataGT,1))';
        toc
        fprintf('save: %s\n',fnameScans);
        save(fnameScans,'dataGT','scan_idxs_gt');
    end
elseif (evalMode == 6 || evalMode == 7)
    fnameScans = '/BS/leonid-pose/work/experiments-caesar/model-mpii-all-samples-warm-sched-3-bootstrap-4/model/nrdAll.mat';
    l = load(fnameScans,'data','scan_idxs');
    dataGT = l.data;
    scan_idxs_gt = l.scan_idxs;
    assert(sum(scan_idxs_gt-scan_idxs) == 0);
end
totalSamples = size(data,1);
dimensions =size(data,2);

arrayScansPointsAux=[];

switch(evalMode)
    case 0
        % Compactness
        [cumulativeStd, evalues, evectors, meanData] = ErrorEvaluation(evalMode,data',arrayScansPointsAux',dimensions,totalSamples,0,totalSamples-1,totalSamples,numComponents);
        saveToDirC = [saveDir '/compactness/'];
        if (~exist(saveToDirC, 'dir'))
            mkdir(saveToDirC);
        end
        save([saveToDirC '/evectors.mat'], 'evectors');
        save([saveToDirC '/meanData.mat'], 'meanData');
        save([saveToDirC '/evalues.mat'], 'evalues');
        save([saveToDirC '/cumulativeStd.mat'], 'cumulativeStd');
        
    case 1
        if (id < totalSamples)
            % Generalization (parallelized): computed distances are saved to calculate mean/std afterwards
%             [distances] = ErrorEvaluation(evalMode,data',arrayScansPointsAux',dimensions,totalSamples,0,totalSamples-1,totalSamples,numComponents,id);
            nPCA = 50;
            if (nTraininMeshes > totalSamples)
                nTraininMeshes = totalSamples;
            end
            
            [distances] = ErrorEvaluation(evalMode,data',arrayScansPointsAux',dimensions,totalSamples,0,nPCA-1,totalSamples,numComponents,id,nTraininMeshes);
            saveToDirG = [saveDir '/generalization/'];
            if (~exist(saveToDirG, 'dir'))
                mkdir(saveToDirG);
            end
            save([saveToDirG '/distances_pcacompidx_' num2str(id) '.mat'], 'distances' );
            fprintf('Num dist: %d\n', size(distances,1));
        end
    case 2
        % Specificity
        sampleSize = 10000;
        nPCA = 50;
        if (nTraininMeshes > totalSamples)
            nTraininMeshes = totalSamples;
        end
%         [meandist, covdist] = ErrorEvaluation(evalMode,data',arrayScansPointsAux',dimensions,totalSamples,0,totalSamples-1,sampleSize,numComponents,id);
        [meandist, covdist] = ErrorEvaluation(evalMode,data',arrayScansPointsAux',dimensions,totalSamples,0,nPCA-1,sampleSize,numComponents,id,nTraininMeshes);
        saveToDirS = [saveDir '/specificity/'];
        if (~exist(saveToDirS, 'dir'))
            mkdir(saveToDirS);
        end
        save([saveToDirS '/covdist' num2str(id) '.mat'], 'covdist');
        save([saveToDirS '/meandist' num2str(id) '.mat'], 'meandist');

    case 3
        % mean/std computation for generalization
        fprintf('Compute mean/cov for generalization\n');
        DistancesTotal=zeros(totalSamples,totalSamples);
        fprintf('Loading distances\n');
        for i=1:totalSamples
            fName = [saveDir '/distances_pcacompidx_' num2str(id) '.mat'];
            if ~exist(fName, 'file')
                fprintf('File missing\n');
                fprintf('%s\n',fName);
                assert(0);
            end
            load(fName);
            DistancesTotal(i,:) = distances;
        end
        	
    meandist = mean(DistancesTotal);
    covdist = diag(cov(DistancesTotal,1));
    %[meandist, covdist] =  ErrorEvaluation(evalMode,DistancesTotal,0,totalSamples,dimensions);
    fprintf('Computing mean/cov\n');
    save([saveDir '/covdist' num2str(id) '.mat'], 'covdist');
    save([saveDir '/meandist' num2str(id) '.mat'], 'meandist');
	
    case 4
        if (id < totalSamples)
            nPCA = 50;
            if (nTraininMeshes > totalSamples)
                nTraininMeshes = totalSamples;
            end
            
%             figure(100); clf;
%             pointsTemplate = reshape(data(1,:)',19347/3,3);
%             pointsScan = reshape(dataGT(1,:)',19347/3,3);
%             plot3(pointsScan(:,1),pointsScan(:,2),pointsScan(:,3),'g+')
%             axis equal; hold on;
%             plot3(pointsTemplate(:,1),pointsTemplate(:,2),pointsTemplate(:,3),'r+')
            
            [distances] = ErrorEvaluation(evalMode,data',arrayScansPointsAux',dimensions,totalSamples,0,nPCA-1,totalSamples,numComponents,id,nTraininMeshes,dataGT');
            saveToDirG = [saveDir '/generalizationGT/'];
            if (~exist(saveToDirG, 'dir'))
                mkdir(saveToDirG);
            end
            save([saveToDirG '/distances_pcacompidx_' num2str(id) '.mat'], 'distances' );
            fprintf('Num dist: %d\n', size(distances,1));
        end
    
    case 5
        % Specificity
        sampleSize = 10000;
        nPCA = 50;
        if (nTraininMeshes > totalSamples)
            nTraininMeshes = totalSamples;
        end
        [meandist, covdist] = ErrorEvaluation(evalMode,data',arrayScansPointsAux',dimensions,totalSamples,0,nPCA-1,sampleSize,numComponents,id,nTraininMeshes,dataGT');
        saveToDirS = [saveDir '/specificityGT/'];
        if (~exist(saveToDirS, 'dir'))
            mkdir(saveToDirS);
        end
        save([saveToDirS '/covdist' num2str(id) '.mat'], 'covdist');
        save([saveToDirS '/meandist' num2str(id) '.mat'], 'meandist');
    
    case 6
        if (id < totalSamples)
            nPCA = 50;
            if (nTraininMeshes > totalSamples)
                nTraininMeshes = totalSamples;
            end

            [distances] = ErrorEvaluation(evalMode-2,data',arrayScansPointsAux',dimensions,totalSamples,0,nPCA-1,totalSamples,numComponents,id,nTraininMeshes,dataGT');
            saveToDirG = [saveDir '/generalizationOurs4/'];
            if (~exist(saveToDirG, 'dir'))
                mkdir(saveToDirG);
            end
            save([saveToDirG '/distances_pcacompidx_' num2str(id) '.mat'], 'distances' );
            fprintf('Num dist: %d\n', size(distances,1));
        end
    
    case 7
        % Specificity
        sampleSize = 10000;
        nPCA = 50;
        if (nTraininMeshes > totalSamples)
            nTraininMeshes = totalSamples;
        end
        [meandist, covdist] = ErrorEvaluation(evalMode-2,data',arrayScansPointsAux',dimensions,totalSamples,0,nPCA-1,sampleSize,numComponents,id,nTraininMeshes,dataGT');
        saveToDirS = [saveDir '/specificityOurs4/'];
        if (~exist(saveToDirS, 'dir'))
            mkdir(saveToDirS);
        end
        save([saveToDirS '/covdist' num2str(id) '.mat'], 'covdist');
        save([saveToDirS '/meandist' num2str(id) '.mat'], 'meandist');    
        
    otherwise
        disp('Nothing to do');
end

fprintf('Done\n');

end