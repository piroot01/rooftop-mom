function G = basisFuncGraph(countX, countY)
%{
    % The graph for the rooftop basis functions on a rectangular grid
    % of size countX x countY patches in each direction will have
    % countX*countY vertices in total.
    vertices = 1:(countX * countY);

    % The graph edges must be labeled in a specifid way; first we will
    % consider the exges (basis functions) that are oriented in the
    % x direction. They have a simple connection list [n, n+1], except
    % the case when mod(n, countX) == 0, then the edge does not exist.
    edgesXFst = vertices;
    edgesXFst(countX:countX:end) = [];
    edgesXSnd = edgesXFst + 1;
    edgesX = [edgesXFst; edgesXSnd];

    % Then we can construct the edges of the graph which are associated
    % with the basis functions in the y direction and they are indexed
    % from the left to right from the bottom to top.
    edgesYFst = vertices(1:(end - countX));
    edgesYSnd = edgesYFst + countX;
    edgesY = [edgesYFst; edgesYSnd];

    edges = [edgesX, edgesY];

    G = graph(edges(1, :), edges(2, :));
%}
    % Faster version of the code above.
    numVertices = countX * countY;

    % Source nodes in the x direction.
    sX = 1:numVertices;
    sX(countX:countX:end) = [];

    % Source nodes in the y direction.
    sY = 1:(numVertices - countX);

    % Concatenate and calculate targets directly inside the constructor
    % sX + 1 gives the x targets, sY + countX gives the y targets
    G = graph([sX, sY], [sX + 1, sY + countX], numVertices);
end
