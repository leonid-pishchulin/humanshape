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