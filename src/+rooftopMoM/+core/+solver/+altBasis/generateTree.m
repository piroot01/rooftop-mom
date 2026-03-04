function tree = generateTree(countX, countY)
% generateTree: Returns a tree basis function as a 1xN matrix.
%
%   Inputs:
%       countX: Positive integer specifying the number of patches in the x direction.
%       countY: Positive integer specifying the number of patches in the y direction.
%
%   Output:
%       tree: 1xN matrix describing the tree.

    % Create the mesh represented by a matrix.
    treeBasis = reshape(1:(countX * countY), countX, countY);

    % Flip the even columns upside down using flipud.
    treeBasis(:, 2:2:end) = flipud(treeBasis(:, 2:2:end));

    % Flatten into a 1xN row vector.
    tree = treeBasis(:).';
end
