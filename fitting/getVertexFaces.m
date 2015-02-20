function vertexFaces = getVertexFaces(vertexIdx, facesAll)

vertexFaces = [];

vertexFaces = vertcat(vertexFaces,find(facesAll(:,1)==vertexIdx));
vertexFaces = vertcat(vertexFaces,find(facesAll(:,2)==vertexIdx));
vertexFaces = vertcat(vertexFaces,find(facesAll(:,3)==vertexIdx));

vertexFaces = vertexFaces';

end