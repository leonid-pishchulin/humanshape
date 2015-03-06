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

function normals = getNormals(vertices, faces)
% Generates the normals for each face given the vertice & face information
% as input.

nFaces = size(faces,1);

normals = zeros(nFaces,3);

i=1:nFaces;
normals(i,:) = cross(vertices(faces(i,1),:)-vertices(faces(i,2),:), vertices(faces(i,1),:)-vertices(faces(i,3),:));
l = sqrt(sum(normals.^2,2));

normals = normals ./ repmat(l,1,3);

end




