function determinant = determinant(obj, content)
% determinant: Returns the determinant of the 2x2 blocks at every frequency, or the single block of a strip.
%
%   Input:
%       content: rooftopMoM.basisFunctions.RooftopContent.
%
%   Output:
%       determinant: Array of the size of the frequency grid.

    if content == rooftopMoM.basisFunctions.RooftopContent.XONLY
        determinant = obj.blockXX_;
    elseif content == rooftopMoM.basisFunctions.RooftopContent.YONLY
        determinant = obj.blockYY_;
    else
        determinant = obj.blockXX_ .* obj.blockYY_ - obj.blockXY_ .* obj.blockYX_;
    end
end
