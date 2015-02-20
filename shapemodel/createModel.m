function createModel(expidx)
fprintf('createModel\n');

p = expParams(expidx);

load([p.modelOutDir '/meanData-flat.mat'],'meanData');
dim = size(meanData,2);
nPoints = dim./3;

% read model mesh
fid = fopen(p.modelMesh,'rt');
A = textscan(fid, '%f', 'HeaderLines',1);
A = A{1};
fclose(fid);
data = reshape(A(1:83837),13,nPoints);
data = data';
% reading 3D points
meanOld = data(:,1:3);
%reading the rest of the line
params = data(:,4:13);
%rest of the values in file
r1 = A(83838:122519);
r2 = A(122520:122719);
data1 = reshape(r1,3,size(r1,1)/3); %8682
data1 = data1';
data2 = reshape(r2,8,size(r2,1)/8);

% mean after PCA
meanNew = reshape(meanData,nPoints,3);

load([p.modelOutDir '/evectors-flat.mat'],'evectors');
totalVectors = size(evectors,2)/dim;
evectors = reshape(evectors,totalVectors,dim);

[regParams,Bfit,ErrorStats]=absor(meanOld',meanNew');

meanNew = meanNew * regParams.R;
meanNew = meanNew + repmat(regParams.t',size(meanOld,1),1);
modelDif = meanOld(1,:) - meanNew(1,:);
gg = repmat(modelDif,6449,1);
meanNew = meanNew+gg;

for i = 1:size(evectors,1)
    ev = reshape(evectors(i,:),size(meanOld,1),size(meanOld,2));
    ev = ev * regParams.R;
    evectors(i,:) = reshape(ev, 1, size(ev,1)*size(ev,2));
end
save([p.modelOutDir '/evectors.mat'],'evectors');

%writting model file with the new aligned mean
fid = fopen([p.modelOutDir '/model.dat'],'w');

%writting the header
fprintf(fid,[num2str(nPoints) ' 12894 25 5 \n']);

j = 1;
while(j <= nPoints)
    finalData(j,1:3) = meanNew(j,:);
    %format shortG
    finalData(j,4:13) = params(j,:);
    fprintf(fid,'%f %f %f %d %f %d %f %d %f %d %d %d %d\n', finalData(j,:));
    j = j + 1;
end
elements=size(data1,1);
j=1;
while(j < elements+1)
    fprintf(fid,'%d %d %d \n', data1(j,:));
    j = j + 1;
end
elements = size(data2,2);
j=1;
while(j < elements + 1)
    fprintf(fid,'%d %f %f %f %f %f %f %d \n', data2(:,j));
    j = j + 1;
end
fclose(fid);

points = meanNew;

% save mean shape
save([p.modelOutDir '/meanShape.mat'],'points');
end