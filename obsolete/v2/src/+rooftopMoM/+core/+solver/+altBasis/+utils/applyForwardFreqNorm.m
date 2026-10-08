function I = applyForwardFreqNorm(I, loopCount, omega)
% applyForwardFreqNorm: Applies a factor of (i * omega) to the currents associated with the tree basis.
%
%   Inputs:
%       I: Nx1 loop-tree current vector.
%       loopCount: Total number of loops on the antenna.
%       omega: Scalar representing the angular frequency of the incident wave.
%
%   Output:
%       I: Nx1 frequency-normalized vector.

    I = rooftopMoM.core.solver.altBasis.utils.applyToTreeCurrent(I, loopCount, @(x) x .* (1i * omega));
end
