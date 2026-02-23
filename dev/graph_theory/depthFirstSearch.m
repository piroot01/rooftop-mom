function [sourceTreeVertices, targetTreeVertices] = depthFirstSearch(G)
    treeVertices = dfsearch(G, 1, 'edgetonew');
    sourceTreeVertices = treeVertices(:, 1);
    targetTreeVertices = treeVertices(:, 2);
end
