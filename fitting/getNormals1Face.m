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