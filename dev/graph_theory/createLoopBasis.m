function loopBasis = createLoopBasis(countX, countY)
% loopBasis: Return 4xN matrix, where each column represents four
%   vertices of a single loop basis function consisting of four
%   rooftop basis functions.
%
%   Inputs:
%       countX: Positive integer specifying the number of patches in the x direction.
%       countY: Positive integer specifying the number of patches in the y direction.
%
%   Output:
%       loopBasis: 4xN matrix describing the loop basis functions.

    % The bottom-left corner of any loop cannot be in the top row.
    % So, we only consider vertices up to the second-to-last row.
    n = 1:(countX * (countY - 1));

    % The bottom-left corner also cannot be on the right-most edge.
    % We remove all indices where mod(n, countX) == 0.
    n(countX:countX:end) = [];

    % Now, 'n' contains all the bottom-left vertices of every single loop.
    % Each loop has its vertices indexed CW.
    loopBasis = [n; n + 1; n + countX + 1; n + countX];
end
