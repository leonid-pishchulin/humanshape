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

function vertexFaces = getVertexFaces(vertexIdx, facesAll)

vertexFaces = [];

vertexFaces = vertcat(vertexFaces,find(facesAll(:,1)==vertexIdx));
vertexFaces = vertcat(vertexFaces,find(facesAll(:,2)==vertexIdx));
vertexFaces = vertcat(vertexFaces,find(facesAll(:,3)==vertexIdx));

vertexFaces = vertexFaces';

end