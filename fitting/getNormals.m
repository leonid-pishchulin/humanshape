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




