function I = applyBackwardFreqNorm(I, loopCount, omega)
% applyBackwardFreqNorm: Applies a factor of (i * omega)^(-1) to the currents associated with the loop basis.
%
%   Inputs:
%       I: Nx1 loop-tree current vector.
%       omega: Scalar representing the angular frequency of the incident wave.
%       loopCount: Total number of loops on the antenna.
%
%   Output:
%       I: Nx1 frequency-normalized vector.

    I = rooftopMoM.core.solver.altBasis.utils.applyToLoopCurrent(I, loopCount, @(x) x ./ (1i * omega));
end
