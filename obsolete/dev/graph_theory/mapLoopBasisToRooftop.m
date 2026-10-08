function matrix = mapLoopBasisToRooftop(loopBasis, countX, countY)
    numBasisX = (countX - 1) * countY;
    offset = numBasisX - countX; 

    % NOTE: The way that the loop basis functions are generated we
    % can avoid using the loopBasis for the x-dir rooftops
    matrix = [
        1:size(loopBasis, 2);
        loopBasis(2, :) + numBasisX;
        countX:numBasisX;
        loopBasis(4, :) + offset
    ];
end
