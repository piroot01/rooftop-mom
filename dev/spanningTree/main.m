%% Fist we generate the primary graph
clc; close all;

countX = 10; % number of patches in the x-direction
countY = 10; % number of patches in the y-direction

% Now let us introduce the indexing scheme:
% For the vertices the scheme will be zig-zag from left bottom corner to the top
% The patches in the same manner.

% To generate the rooftops, we can just consider the dual mesh which is the same
% excluding the boundary edges, therefore its dimensions are (countX - 1) x (countY - 1)
% We must point out that the indices are now patch indices of the primary mesh, but vertices
% of the dual mesh.
[~, rooftops] = generateRectMeshGraph(countX - 1, countY - 1);

% The dual basis made of BC (Buffa-Christiansen) functions live on the dual grid, but the patches
% are vertices of the primary mesh without the boundary vertices (all vertices that coincide with
% a boundary edge).
%[bCVertices, buffaChristiansen, info] = generateRectMeshGraph(countX, countY, true);

% Main grid
[mainVerts, mainEdges, mainInfo] = generateRectMeshGraph(countX, countY);

% Extract the internal sub-grid
[interVerts, interEdges, interInfo] = extractInteriorGrid(mainVerts, mainEdges, countX + 1, countY + 1);

% Choose random point on the internal mesh boundary.
interTreeSeedVertex = interInfo.boundaryVertices(randi([1, numel(interInfo.boundaryVertices)], 1));

% Generate full mesh adjacency matrix.
mainAdjMatrix = generateAdjMatrix(mainVerts, mainEdges);

% Generate the adjacency matrix of the internal mesh.
mainVertCount = numel(mainVerts);
interAdjMatrix = sparse(mainVertCount, mainVertCount);
% Reindex
interAdjMatrix(interVerts, interVerts) = mainAdjMatrix(interVerts, interVerts);

% Generate the interior spanning tree.
interTree = generateSpanningTreeDFS(interAdjMatrix, interTreeSeedVertex);

% Crop the zero elements of the interior tree.
interVertCount = numel(interVerts);
interTree = interTree(1:interVertCount - 1, :);

% Prepend the tree witth the intial internal tree edge.
firstInterVert = intersect(mainVerts(mainAdjMatrix(:, interTreeSeedVertex)), mainInfo.boundaryVertices);
% If more neighboring vertices are available, choose random one.
firstInterVert = firstInterVert(randi([1, numel(firstInterVert)], 1));
interTree = cat(1, [firstInterVert, interTree(1, 1)], interTree);

% We need to implement following functions:
% 1. Returns a list of vertices that are coindicent of a given vertex
% 2. Returns true if given rooftop (edge) crosses BC basis
% 3. Verifies if the given rooftop-tree is a spanning tree.
% 4. Verifies if the given BC-tree is a spanning tree.
% 5. Verifies if given rooftop-tree and BC-tree forms a complementary pair.
%   - the rooftop-tree must have every vertex of its graph explored
%   - the BC-tree mush have every internal vertex of the primary mesh graph
%       explored
%   - there must not be crossing pair of rooftop and BC function

%{
% Simple pipeline for a tree generation
tic;
[vert, edg] = generateRectMeshGraph(countX, countY);
toc;

tic;
adjMat = generateAdjMatrix(vert, edg);
toc;

tic;
tree = generateSpanningTreeDFS(adjMat, randi([1, numel(vert)], 1));
toc;
%}

% Function that extracts the interior sub-graph from an M x N grid.
function [intV, intE, info] = extractInteriorGrid(V, E, countX, countY)
    % Map 1D vertex indices to 2D coordinates (xs, ys).
    xs = mod(V - 1, countX) + 1;
    ys = floor((V - 1) / countX) + 1;

    % Find vertices that are strictly inside the boundary
    % They must be greater than 1, and less than the max dimensions.
    isInterior = (xs > 1) & (xs < countX) & (ys > 1) & (ys < countY);

    % Extract the interior vertices (retaining their original graph
    % indices).
    intV = V(isInterior);

    % xtract the boundary of the internal grid:
    % The internal boundary is the outermost ring of the interior vertices.
    isInternalBoundary = isInterior & ...
        (xs == 2 | xs == (countX - 1) | ys == 2 | ys == (countY - 1));
    info.boundaryVertices = V(isInternalBoundary);

    % Filter the edges:
    % Create a boolean lookup array for the interior status of any node.
    isInteriorLookup = false(max(V), 1);
    isInteriorLookup(intV) = true;

    validEdges = isInteriorLookup(E(:, 1)) & isInteriorLookup(E(:, 2));
    intE = E(validEdges, :);
end

% Generates a symmetric sparse adjacency matrix for an undirected graph
function adjMatrix = generateAdjMatrix(V, E)
    vertexCount = numel(V);
    
    % Concatenate both directions to ensure the matrix is symmetric
    fromNodes = [E(:, 1); E(:, 2)];
    toNodes   = [E(:, 2); E(:, 1)];
    
    adjMatrix = sparse(fromNodes, toNodes, true, vertexCount, vertexCount);
end

% Function that for given graph performs a iterative DFS spanning tree generation
function treeEdges = generateSpanningTreeDFS(adjMatrix, seedVertex)
    vertexCount = size(adjMatrix, 1);
    visitedVertices = false(1, vertexCount);

    % Preallocate the stack.
    stack = zeros(vertexCount, 1);
    stackSize = 1;
    stack(1) = seedVertex;

    % Mark the seed vertex as visited.
    visitedVertices(seedVertex) = true;

    % Preallocate the tree.
    treeEdges = zeros(vertexCount - 1, 2);
    edgeCount = 0;

    while stackSize > 0
        % Obtain the vertex from the top of the stack.
        currentVertex = stack(stackSize);

        adjVertices = find(adjMatrix(:, currentVertex));

        % Filter out the ones which are explored.
        unvisitedAdjVertices = adjVertices(~visitedVertices(adjVertices));

        if isempty(unvisitedAdjVertices)
            % Dead branch, so we decrement the stackSize.
            stackSize = stackSize - 1;
        else
            % Branch: we choose a random unvisited adjacent vertex.
            nextVertex = unvisitedAdjVertices(randi([1, numel(unvisitedAdjVertices)], 1));
            %nextVertex = max(unvisitedAdjVertices);
            %nextVertex = min(unvisitedAdjVertices);

            % Mark the nextVertex as visited.
            visitedVertices(nextVertex) = true;

            % Save the nextVertex as the next tree edge.
            edgeCount = edgeCount + 1;
            treeEdges(edgeCount, :) = [currentVertex, nextVertex];

            % Increment the stack.
            stackSize = stackSize + 1;
            stack(stackSize) = nextVertex;
        end
    end
end

% Function that returns a graph of a rectangular 2D mesh grid of size
% countX by countY patches.
function [V, E, info] = generateRectMeshGraph(countX, countY, removeBoundaryEdges)
    info = [];

    if nargin < 3
        removeBoundaryEdges = false;
    end

    % Calculate the total number of vertices in each direction.
    vertexCountX = countX + 1;
    vertexCountY = countY + 1;
    vertexCount = vertexCountX * vertexCountY;

    V = (1:vertexCount).';

    % Edges in the x-direction (horizontal):
    % Connect nodes to their rightward neighbor (k, k + 1).
    % We exclude the right-most node of each row using modulo arithmetic 
    % so edges don't incorrectly wrap around to the start of the next row.
    startNodesX = V(mod(V, vertexCountX) ~= 0);

    if removeBoundaryEdges
        % Remove top boundary (first row) and bottom boundary (last row)
        % Keep only start nodes that are strictly greater than the first row 
        % and less than or equal to the start of the last row.
        startNodesX = startNodesX(startNodesX > vertexCountX & startNodesX <= vertexCountX * countY);

        info.internalVertices = startNodesX(2:end);
    end

    edgeX = [startNodesX, startNodesX + 1];

    % Edges in the y-direction (vertical):
    % Connect nodes to their downward/upward neighbor (k, k + vertexCountX).
    % We use all nodes except those in the very last row as starting points.
    startNodesY = (1:(vertexCountX * countY))';

    if removeBoundaryEdges
        % Remove left boundary (first column) and right boundary (last column)
        % Keep only start nodes where the modulo is neither 1 (left) nor 0 (right).
        startNodesY = startNodesY(mod(startNodesY, vertexCountX) ~= 1 & mod(startNodesY, vertexCountX) ~= 0);
    end

    % Store the boundary vertices.
    boundaryVerticesY = (vertexCountX + 1):vertexCountX:(vertexCount - vertexCountX);
    info.boundaryVertices = [1:vertexCountX, (vertexCount - vertexCountX + 1):vertexCount, ...
        boundaryVerticesY, ...
        boundaryVerticesY + vertexCountX - 1].';

    edgeY = [startNodesY, startNodesY + vertexCountX];

    E = [edgeX; edgeY];
end


%% Plot the grid
sizeX = 1; % size of the single patch in the x-direction
sizeY = 1; % and the in the y-direction

% Generate the coordinates
[xs, ys] = meshgrid(0:sizeX:(countX * sizeX), 0:sizeY:(countY * sizeY));

% CRITICAL FIX: Transpose the matrices so that 1D indexing (xs(:)) 
% moves left-to-right first, matching your edge generation logic!
xs = xs';
ys = ys';

% Now we would like to make a dual graph mesh that will be for the rooftops.
[xsRoofops, ysRooftop] = meshgrid((0.5 * sizeX):sizeX:(countX * sizeX), ...
    (0.5 * sizeY):sizeY:(countY * sizeY));

xsRoofops = xsRoofops';
ysRooftop = ysRooftop';

figure;
% Plot the grid vertices
plot(xs(:), ys(:), "g.", "MarkerSize", 10, "LineWidth", 2);
hold on;

% Plot the rooftop (dual) vertices
plot(xsRoofops(:), ysRooftop(:), "y.", "MarkerSize", 10, "LineWidth", 2);

% Get the 1D indices for the start and end of every edge
startNodes = interTree(:, 1);
endNodes = interTree(:, 2);

% Extract the actual X and Y coordinates for those nodes
xStarts = xs(startNodes);
xEnds = xs(endNodes);
yStarts = ys(startNodes);
yEnds = ys(endNodes);

% When plotted, MATLAB stops drawing at a NaN, allowing us to 
% draw thousands of separate line segments in one incredibly fast command.
xNan = [xStarts, xEnds, nan(size(xStarts))]';
yNan = [yStarts, yEnds, nan(size(yStarts))]';

plot(xNan(:), yNan(:), 'w-', 'LineWidth', 3);

axis equal;
hold off;
