classdef Constants
% Constants: Physical constants of the vacuum in SI units (CODATA 2018).
%   The speed of light c and the impedance Z0 are given; mu0 = Z0 / c and epsilon0 = 1 / (Z0 * c) are derived from
%   them, so c^2 * mu0 * epsilon0 = 1 and Z0 = sqrt(mu0 / epsilon0) hold up to rounding.

    properties (Constant)
        % Speed of light in vacuum (m/s), exact.
        c = 299792458;

        % Impedance of free space (ohm).
        Z0 = 376.730313668;

        % Magnetic constant (H/m).
        mu0 = rooftopMoM.assembler.Constants.Z0 / rooftopMoM.assembler.Constants.c;

        % Electric constant (F/m).
        epsilon0 = 1 / (rooftopMoM.assembler.Constants.Z0 * rooftopMoM.assembler.Constants.c);
    end
end
