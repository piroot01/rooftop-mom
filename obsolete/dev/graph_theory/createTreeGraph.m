function treeG = createTreeGraph(G, sourceTreeVertices, targetTreeVertices)
    treeG = graph(sourceTreeVertices, targetTreeVertices, G.numnodes);
end
