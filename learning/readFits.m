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

function data = readFits(fitDir)

fprintf('readFits()\n');

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
