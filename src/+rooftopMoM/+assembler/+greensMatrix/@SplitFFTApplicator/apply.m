function output = apply(obj, input)
% apply: Returns the product of the Green's matrix with a vector of rooftop coefficients.
%
%   Input:
%       input: Nx1 vector of the rooftop coefficients.
%
%   Output:
%       output: Nx1 vector of the product.

    switch obj.content_
        case rooftopMoM.basisFunctions.RooftopContent.BOTH
            output = obj.applyBothDirections_(input);
        case rooftopMoM.basisFunctions.RooftopContent.XONLY
            output = obj.applySingleDirection_(input, obj.spectra_(1).blockXX(), obj.spectra_(2).blockXX(), ...
                obj.modulationX_, obj.demodulationX_);
        otherwise
            output = obj.applySingleDirection_(input, obj.spectra_(1).blockYY(), obj.spectra_(2).blockYY(), ...
                obj.modulationY_, obj.demodulationY_);
    end
end
