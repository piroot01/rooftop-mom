function loopGraph = createLoopGraph(G, treeG)
    edgesG = G.Edges.EndNodes;
    edgesTreeG = treeG.Edges.EndNodes;

    edgesGSorted = sort(edgesG, 2);
    edgesTreeGSorted = sort(edgesTreeG, 2);

    [~, nonTreeVertexIdx] = setdiff(edgesGSorted, edgesTreeGSorted, 'rows');

    chords = edgesG(nonTreeVertexIdx, :);

    loopGraph = graph(chords(:, 1), chords(:, 2), G.numnodes);
end
