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

function NN = getVertexNN(faces)

[countFaces,~] = size(faces);
countVertices = max(max(faces));

NN = cell(countVertices,1);

for i = 1:countFaces
    
    NN{faces(i,1),1} = horzcat(NN{faces(i,1),1}, [faces(i,2) faces(i,3)]);
    NN{faces(i,2),1} = horzcat(NN{faces(i,2),1}, [faces(i,1) faces(i,3)]);
    NN{faces(i,3),1} = horzcat(NN{faces(i,3),1}, [faces(i,1) faces(i,2)]);
    
end

for i = 1:countVertices
    
    NN{i,1} = unique(NN{i,1});
    NN{i,1} = sort(NN{i,1},'descend');
    
end

end