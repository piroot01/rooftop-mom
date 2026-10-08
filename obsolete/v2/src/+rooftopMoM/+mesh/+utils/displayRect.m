function displayRect(vertices, faces, genome)
% displayMesh: Display a 2D quad mesh and color the active faces.
%
%   Inputs:
%       vertices: 3×N array, each column is a 3D vertex coordinate.
%       faces: 4×F integer array, each column lists vertex indices of a quad face.
%       genome: (optional) Tells which faces to color.

    faceCount = size(faces, 2);

    % Default to the first face if unspecified.
    if nargin < 3 || isempty(genome)
        genome = true(faceCount, 1);
    end

    % Project the vertices to 
    [u, v] = rooftopMoM.mesh.utils.generateRectLocalCoords(vertices, faces);

    localVertices = rooftopMoM.mesh.utils.generateRectLocalVertices(vertices, u, v);

    figure('Name', 'displayRect');
    hold on;

    % Loop through each face and draw it.
    for n = 1:faceCount
        F = faces(:, n);
        V = localVertices(:, F);

        % Choose color based on genome
        if genome(n)
            % black
            fc = [0, 0, 0];
        else
            % white
            fc = [1, 1, 1];
        end

        patch('XData', V(1, :), 'YData', V(2, :), ...
              'FaceColor', fc, 'EdgeColor', [0, 0, 0]);
    end

    axis equal;
    axis off;
    hold off;
end
