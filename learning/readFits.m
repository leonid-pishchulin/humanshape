function data = readFits(fitDir)

fprintf('readNRD()\n');

files = dir(fitDir);
bIsFitDir = zeros(length(files),1);
for i=1:length(files)
    bIsFitDir(i) = exist([fitDir '/' files(i).name '/NRD.mat'],'file');
end
files(bIsFitDir == 0) = [];

load([fitDir '/' files(1).name '/NRD.mat'],'points');
dim = size(points,1)*size(points,2);
data = zeros(length(files),dim);

scannum = 0;
for i = 1:length(files)
    fprintf('.');
    fname = [fitDir '/' files(i).name  '/NRD.mat'];
    assert(exist(fname, 'file') > 0);
    scannum = scannum + 1;
    load(fname, 'points');
    data(scannum,:) = reshape(points, size(points,1)*size(points,2),1);
    if (~mod(i, 100))
        fprintf(' %d/%d\n',i,length(files));
    end
end
assert(scannum == length(files));
fprintf(' done\n');
end
