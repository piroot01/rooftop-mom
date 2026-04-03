function [transformation, data] = generateLoopTreeToRooftopTransformation(countX, countY, sizeX, sizeY)
% generateRooftopToLoopTreeTransformation: Assembles the sparse transformation matrix that converts
%   currents represented in the loop-tree basis into currents in the rooftop basis.
%
%   Inputs:
%       countX, countY: Grid dimensions of the antenna topology.
%       sizeX, sizeY: Single patch size.
%
%   Output:
%       transformation: NxN sparse matrix, where N is the total number of rooftops.
%       data: Struct that contains export variables for further use.

    rooftopCount = rooftopMoM.basisFunctions.utils.generateRooftopCount(countX, countY);

    % Loop basis
    loops = rooftopMoM.core.solver.altBasis.generateLoops(countX, countY);
    loopsConverted = rooftopMoM.core.solver.altBasis.convertLoopToRooftopIndices(loops, countX, countY);
    loopTransformation = rooftopMoM.core.solver.altBasis.generateRooftopToLoopTransformation(loopsConverted, sizeX, sizeY);

    % Tree basis
    tree = rooftopMoM.core.solver.altBasis.generateTree(countX, countY);
    treeConverted = rooftopMoM.core.solver.altBasis.convertTreeToRooftopIndices(tree, countX, countY);
    treeTransformation = rooftopMoM.core.solver.altBasis.generateRooftopToTreeTransformation(treeConverted, rooftopCount);

    % Transpose and combine to transform the coordinates (currents), not the basis itself.
    transformation = cat(2, loopTransformation.', treeTransformation.');

    % Save the tree into the data.
    data.tree = tree;
    data.loopCount = size(loops, 2);
end
