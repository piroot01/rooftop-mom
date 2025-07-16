classdef Constants
    properties (Constant)
        % Speed of light in vacuum (m/s)
        c = 299792458;

        % Magnetic constant (absolute permeability of vacuum) (H/m)
        mu0 = 4 * pi * 1e-7;

        % Electric constant (absolute permittivity of vacuum) (F/m)
        % CODATA 2018: 8.8541878128e-12 F/m
        epsilon0 = 8.8541878128e-12;

        % Impedance of free space (ohms)
        % CODATA 2018: 376.730313668 ohms
        Z0 = 376.730313668;
    end
end
