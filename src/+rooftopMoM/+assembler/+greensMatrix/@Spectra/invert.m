function inverse = invert(obj, content)
% invert: Returns the spectra of the inverse of the block circulant matrix. The FFT diagonalizes every circulant
%   block, so the inverse is the inverse of the 2x2 matrix [XX, XY; YX, YY] at every frequency, the reciprocal of the
%   single block for a strip.
%
%   Input:
%       content: rooftopMoM.basisFunctions.RooftopContent.
%
%   Output:
%       inverse: rooftopMoM.assembler.greensMatrix.Spectra.

    if content == rooftopMoM.basisFunctions.RooftopContent.XONLY
        inverse = rooftopMoM.assembler.greensMatrix.Spectra(1 ./ obj.blockXX_, [], [], []);
    elseif content == rooftopMoM.basisFunctions.RooftopContent.YONLY
        inverse = rooftopMoM.assembler.greensMatrix.Spectra([], [], [], 1 ./ obj.blockYY_);
    else
        determinant = obj.determinant(content);
        inverse = rooftopMoM.assembler.greensMatrix.Spectra(obj.blockYY_ ./ determinant, ...
            -obj.blockXY_ ./ determinant, -obj.blockYX_ ./ determinant, obj.blockXX_ ./ determinant);
    end
end
