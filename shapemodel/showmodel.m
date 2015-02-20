%
% showmodel(vertices, faces, color, [facelist], [surface], [colors])
%
% faces - for every face coordinates of the corners as column vectors
%             [v1 v2 v3]
%
function vertices = showmodel(vertices, faces, color, facelist, surface, colors)

if ~exist('surface','var')
	surface = 2;
end

if ~exist('color','var')
    color = 'b';
end

if size(vertices,1) < size(vertices,2)
    vertices = vertices';
end

% assume that vertices consists of blocks of 3 each describing one triangle
if ~exist('faces','var') || isempty(faces)
    faces = reshape(1:size(vertices,1), 3,size(vertices,1)/3)';
    for i = 1:3:size(vertices,1)
        vertices(i:i+2,:) = vertices(i:i+2,:)';
    end
end

if surface == 1
	vertices(:,3) = vertices(:,3) - min(vertices(:,3));
elseif surface == 2
	vertices = upright(vertices);
	vertices(:,3) = vertices(:,3) - min(vertices(:,3));
elseif surface == 3
	vertices = upright(vertices, [83, 4242, 985, 3300, 1116]);
	vertices(:,3) = vertices(:,3) - min(vertices(:,3));
elseif surface == 4
    % do nothing
elseif length(surface) > 1
	vertices = upright(vertices, surface);
	vertices(:,3) = vertices(:,3) - min(vertices(:,3));
elseif surface < 0
	vertices = -upright(vertices);
	vertices(:,3) = vertices(:,3) - min(vertices(:,3));
end

if exist('facelist','var') && ~isempty(facelist)
    f = faces(facelist,:);
else
    f = faces;
end

x = reshape(vertices(f',1), 3, size(f,1));
% matlab's coordinate system is lefthanded
y = reshape(vertices(f',2), 3, size(f,1));
z = reshape(vertices(f',3), 3, size(f,1));

c = mean(vertices);
%axis([c(1)-1 c(1)+1 c(2)-1 c(2)+1 c(3)-1 c(3)+1])
axis on
grid on
%shading interp
h=light;
lightangle(h,0,30);
h=light;
lightangle(h,120,30);
h=light;
lightangle(h,240,30);
set(gcf,'renderer','OpenGL');


if exist('colors','var')
    colormap jet;
    patch('Vertices',vertices, 'Faces',f, 'FaceVertexCData',colors, 'FaceColor','flat', 'FaceLighting','gouraud');
    lighting none;
else
    patch(x, y, z, color);
    lighting gouraud
    material shiny
end
