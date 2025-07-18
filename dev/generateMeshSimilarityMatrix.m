clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

% First mesh
countX = 3;
countY = 3;
sizeX = 1;
sizeY = 1;

rectPrimitive1 = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect1 = rooftopMoM.mesh.UniformRect(rectPrimitive1);
mesh1 = rect1.generateMesh();

% Second mesh
countX = 20;
countY = 20;
sizeX = 1;
sizeY = 1;

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

% Obtain the 3x(4*N) array
function rects = filterRects(vertices, faces, faceIndices)
    % NOTE: the vertices can be 2D or 3D it does not matter
    %rects = arrayfun(@(idx) vertices(:, faces(:, idx)), faceIndices, 'UniformOutput', false);
    vertexIndices = faces(:, faceIndices(:));
    rects = vertices(:, vertexIndices(:));
end

% TODO: The function generateRectSize uses the same internal commands so we extract them into separate function
function points = projectVerticesToLocalCoords(vertices, faces, faceIndex)
    if nargin < 3 || isempty(faceIndex)
        faceIndex = 1;
    end

    % Compute local basis vectors.
    [u, v] = rooftopMoM.mesh.utils.generateRectLocalCoords(vertices, faces, faceIndex);

    points = zeros(2, size(vertices, 2));
    % Project vertices onto local axes.
    points(1, :) = u.' * vertices;
    points(2, :) = v.' * vertices;
end

function area = generateRectIntersectionArea(faceVertices)
    %RECTINTERSECTIONAREA  Area of the common intersection of N axis‑aligned rects
    %
    %   area = rectIntersectionArea(FV) returns the area of the
    %   intersection of N rectangles whose vertices are packed
    %   column‑wise in FV (3×(4·N)):
    %      FV(:,1 + 4*(i-1)) = [u_min; v_min; _]   (bottom‑left)
    %      FV(:,2 + 4*(i-1)) = [u_max; v_min; _]   (bottom‑right)
    %      FV(:,3 + 4*(i-1)) = [u_max; v_max; _]   (top‑right)
    %      FV(:,4 + 4*(i-1)) = [u_min; v_max; _]   (top‑left)
    %
    %   Only the first two rows (U, V) are used; the third row
    %   is ignored.  If there is no overlap, returns 0.
    %

    % NOTE: The indices of the faces must be coplanar
    % grab all U and V coords:
    U = faceVertices(1, :);
    V = faceVertices(2, :);

    % each rect contributes 4 cols: [BL, BR, TR, TL]
    lefts   = U(1:4:end);  % bottom‑left  u
    rights  = U(2:4:end);  % bottom‑right u
    bottoms = V(1:4:end);  % bottom‑left  v
    tops    = V(4:4:end);  % top‑left     v

    % compute overlap lengths
    w = min(rights) - max(lefts);
    h = min(tops)   - max(bottoms);

    % if either is negative, no overlap
    area = max(w, 0) * max(h, 0);
end

function matrix = generateMeshSimilarityMatrixTMP(masterMesh, inputMesh, meshPair)
    if nargin < 3 || isempty(meshPair)
        meshPair = [1, 1];
    end

    % First we select the meshes.
    selectedMasterMesh = masterMesh.get(meshPair(1));
    selectedInputMesh = inputMesh.get(meshPair(2));

    masterVertices = selectedMasterMesh.vertices{1};
    masterFaces = selectedMasterMesh.faces{1};

    inputVertices = selectedInputMesh.vertices{1};
    inputFaces = selectedInputMesh.faces{1};

    masterVerticesLocal = projectVerticesToLocalCoords(masterVertices, masterFaces);
    inputVerticesLocal = projectVerticesToLocalCoords(inputVertices, inputFaces);

    countMaster = generateRectFaceCount(masterVertices, masterFaces);
    countInput = generateRectFaceCount(inputVertices, inputFaces);

    matrix = zeros(countInput, countMaster);

    for m = 1:countMaster
        for n = 1:countInput
            rects = [filterRects(masterVerticesLocal, masterFaces, m), filterRects(inputVerticesLocal, inputFaces, n)];
            matrix(n, m) = generateRectIntersectionArea(rects);
        end
    end
end

function count = generateRectFaceCount(vertices, faces)
    [c1, c2] = rooftopMoM.mesh.utils.generateRectGridSize(vertices, faces);
    count = c1 * c2;
end

meshAligned = generateAlignedMesh(mesh1, mesh2);
matrix = generateMeshSimilarityMatrixTMP(mesh1, meshAligned);

prevGenome = ones(9, 1);
prevGenome([2, 4, 6, 8]) = 0;

newGenome = matrix * prevGenome;

save('genome.mat', 'newGenome');
