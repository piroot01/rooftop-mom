function loopSubgraphs = findFourVertexLoops(G)
    % Extracts all 4-vertex cyclic subgraphs (squares) from a mesh graph
    
    A = adjacency(G);
    
    % A^2 gives the number of length-2 paths between any two nodes.
    A2 = A^2;
    
    % Find pairs of vertices (u, v) that share 2 or more common neighbors.
    % We use triu(..., 1) to only look at the upper triangle so we don't 
    % check the same pair twice (i.e., we check u,v but skip v,u).
    [u_nodes, v_nodes] = find(triu(A2, 1) >= 2);
    
    rawLoops = [];
    
    for k = 1:length(u_nodes)
        u = u_nodes(k);
        v = v_nodes(k);
        
        % Find the exact common neighbors using a logical AND on their adjacency rows
        common = find(A(:, u) & A(:, v));
        
        % If they share at least 2 neighbors, we found a 4-cycle!
        if length(common) >= 2
            % If it's a standard grid, length will be exactly 2. 
            % We take the combinations just in case it's a more complex 3D mesh.
            combos = nchoosek(common, 2);
            for c = 1:size(combos, 1)
                % Store the 4 vertices that make up this loop
                rawLoops(end+1, :) = [u, combos(c,1), v, combos(c,2)]; %#ok<AGROW>
            end
        end
    end
    
    % A square has two diagonals. Therefore, the same 4 vertices will be 
    % discovered twice (once for the top-left/bottom-right diagonal pair, 
    % and once for the top-right/bottom-left diagonal pair).
    
    % We sort the vertex indices in each row so they are uniform, 
    % then use 'unique' to filter out the duplicate discoveries.
    if ~isempty(rawLoops)
        loopSubgraphs = unique(sort(rawLoops, 2), 'rows');
    else
        loopSubgraphs = [];
    end
end
