clc;

%% Fist we generate the primary graph
clc;

countX = 4; % number of patches in the x-direction
countY = 2; % number of patches in the y-direction

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
[~, allBuffaChristiansen] = generateRectMeshGraph(countX, countY, true)



% Function that returns a graph of a rectangular 2D mesh grid.
% V: Vector of vertex indices (1 to total vertices).
% E: Nx2 matrix of edges connecting the vertices.
function [V, E] = generateRectMeshGraph(countX, countY, removeBoundaryEdges)
    if nargin < 3
        removeBoundaryEdges = false;
    end

    % Calculate the total number of vertices in each direction.
    vertexCountX = countX + 1;
    vertexCountY = countY + 1;
    totalVertices = vertexCountX * vertexCountY;

    V = (1:totalVertices).';

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

    edgeY = [startNodesY, startNodesY + vertexCountX];

    E = [edgeX; edgeY];
end


%% Plot the grid
sizeX = 1; % size of the single patch in the x-direction
sizeY = 1; % and the in the y-direction

[xs, ys] = meshgrid(0:sizeX:(countX * sizeX), 0:sizeY:(countY * sizeY));

% Now we would like to make a dual graph mesh that will be for the
% rooftops.
[xsRoofops, ysRooftop] = meshgrid((0.5 * sizeX):sizeX:(countX * sizeX), ...
    (0.5 * sizeY):sizeY:(countY * sizeY));

figure;
plot(xs(:), ys(:), "g*", "MarkerSize", 10, "LineWidth", 2);
hold on;
plot(xsRoofops(:), ysRooftop(:), "y*", "MarkerSize", 10, "LineWidth", 2);
axis equal;
