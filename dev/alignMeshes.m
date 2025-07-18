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

% Second mesh
countX = 4;
countY = 6;
sizeX = 1.5;
sizeY = 4;

rectPrimitive2 = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect2 = rooftopMoM.mesh.UniformRect(rectPrimitive2);
mesh2 = rect2.generateMesh();

% NOTE: We will consider only those meshes that are coplanar and share the same origin [0; 0; 0]
function mesh = generateAlignedMesh(masterMesh, inputMesh, meshPair)
    if nargin < 3 || isempty(meshPair)
        meshPair = [1, 1];
    end

    selectedMesh = inputMesh.get(meshPair(2));

    % The masterMesh will not be resized, used just for size extraction
    [c1, c2] = rooftopMoM.mesh.utils.generateRectGridSize(selectedMesh.vertices{1}, selectedMesh.faces{1});

    mesh = generateAlignedMeshCounts(masterMesh, c1, c2, meshPair(1));
end

function mesh = generateAlignedMeshCounts(masterMesh, c1, c2, meshIndex)
    if nargin < 4 || isempty(meshIndex)
        meshIndex = 1;
    end

    selectedMesh = masterMesh.get(meshIndex);

    [s1, s2] = rooftopMoM.mesh.utils.generateRectSize(selectedMesh.vertices{1}, selectedMesh.faces{1});

    newSizeX = s1 / double(c1)
    newSizeY = s2 / double(c2)
    rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(c1, c2, newSizeX, newSizeY);
    rect = rooftopMoM.mesh.UniformRect(rectPrimitive);
    mesh = rect.generateMesh();
end

meshAligned = generateAlignedMesh(mesh1, mesh2);

[s1, s2] = rooftopMoM.mesh.utils.generateRectSize(mesh1.vertices{1}, mesh1.faces{1})
[s1, s2] = rooftopMoM.mesh.utils.generateRectSize(meshAligned.vertices{1}, meshAligned.faces{1})
[c1, c2] = rooftopMoM.mesh.utils.generateRectGridSize(meshAligned.vertices{1}, meshAligned.faces{1})
