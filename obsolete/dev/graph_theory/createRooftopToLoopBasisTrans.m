function matrix = createRooftopToLoopBasisTrans(loopBasis, sizeX, sizeY)
    % The number of the rooftop basis functions will always be
    % contained inside the loopBasis matrix at the second row
    % last column
    numRooftops = loopBasis(2, end);
    numLoops = size(loopBasis, 2);
    %{
    matrix = zeros(numLoops, numRooftops);
    for m = 1:numLoops
        for n = 1:numRooftops
            if loopBasis(1, m) == n
                el = 1 / sizeY;
            elseif loopBasis(2, m) == n
                el = 1 / sizeX;
            elseif loopBasis(3, m) == n
                el = (-1) / sizeY;
            elseif loopBasis(4, m) == n
                el = (-1) / sizeX;
            else
                el = 0;
            end

            matrix(m, n) = el;
        end
    end
    %}

    % optimised variant using sparse matrix
    % Duplicate row coordinates for all 4 edges
    rowIdx = repmat(1:numLoops, 4, 1);
    colIdx = loopBasis(1:4, :);

    % Expand the values to match the size of our coordinates
    vals = repmat([1/sizeY; 1/sizeX; -1/sizeY; -1/sizeX], 1, numLoops);

    % Filter out invalid/zero indices (e.g., boundary edges)
    valid = colIdx > 0 & colIdx <= numRooftops;

    % Assemble the sparse matrix in one shot
    matrix = sparse(rowIdx(valid), colIdx(valid), vals(valid), numLoops, numRooftops);
end
