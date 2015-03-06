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

function normals1Face = getNormals1Face(vertexIdxs,faces,normals)

if (isempty(faces))
    normals1Face = zeros(length(vertexIdxs),0);
    return;
end

nPoints = length(vertexIdxs);
normals1Face = zeros(nPoints,3);
for i = 1:nPoints
    facesPoint = getVertexFaces(vertexIdxs(i), faces);
    if (isempty(facesPoint))
        normal = [NaN NaN NaN];
    else
        normal = normals(facesPoint(1),:);
    end
    normals1Face(i,:) = normal; % picking only the first face
end

end