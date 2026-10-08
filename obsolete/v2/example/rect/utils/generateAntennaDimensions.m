function [countX, countY, sizeX, sizeY] = generateAntennaDimensions(electricalSize, sideRatio, rooftopCountPerOneWavelength)
    if nargin < 3 || isempty(rooftopCountPerOneWavelength)
        rooftopCountPerOneWavelength = 10;
    end

    % We will assume that sX = sY = 1 for now
    sizeX = 1;
    sizeY = 1;

    sideRatio = floor(sideRatio);

    % Then we will set the sorter side to be the one in the x direction
    countX = electricalSize * rooftopCountPerOneWavelength / (pi * sqrt(1 + sideRatio^2));
    countX = ceil(countX);

    countY = sideRatio * countX;
end
