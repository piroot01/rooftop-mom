function output = padWithZeros(input, sizeX, sizeY)
% padWithZeros: Returns the array zero-padded (or cut) to sizeX x sizeY, starting at the top-left corner.
%
%   Inputs:
%       input: 2D array (may be empty).
%       sizeX: Number of rows of the output.
%       sizeY: Number of columns of the output.
%
%   Output:
%       output: sizeX x sizeY array of the class of the input, empty for an empty input.

    if isempty(input)
        output = [];
        return;
    end

    output = zeros(sizeX, sizeY, 'like', input);
    rows = 1:min(size(input, 1), sizeX);
    columns = 1:min(size(input, 2), sizeY);
    output(rows, columns) = input(rows, columns);
end
