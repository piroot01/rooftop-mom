clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

% First mesh
countX = 3;
countY = 3;
sizeX = 2;
sizeY = 1;

rectPrimitive1 = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect1 = rooftopMoM.mesh.UniformRect(rectPrimitive1);
mesh1 = rect1.generateMesh();

function rects = filterRectsTMP(vertices, faces, faceIndices)
    %rects = arrayfun(@(idx) vertices(:, faces(:, idx)), faceIndices, 'UniformOutput', false);
    vertexIndices = faces(:, faceIndices(:));
    rects = vertices(:, vertexIndices(:));
end

indices = [1, 2];

rects = filterRectsTMP(mesh1.vertices{1}, mesh1.faces{1}, indices)
