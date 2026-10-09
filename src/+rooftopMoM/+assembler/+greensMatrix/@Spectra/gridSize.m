function [sizeX, sizeY] = gridSize(obj)
% gridSize: Returns the size of the frequency grid.

    blocks = {obj.blockXX_, obj.blockYY_, obj.blockXY_, obj.blockYX_};
    [sizeX, sizeY] = size(blocks{find(~cellfun(@isempty, blocks), 1)});
end
