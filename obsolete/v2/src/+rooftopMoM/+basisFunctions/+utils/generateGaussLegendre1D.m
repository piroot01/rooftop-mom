function [points, weights] = generateGaussLegendre1D(order, a, b)
% generateGaussLegendre1D: Simple wrapper for Greg von Winckel’s lgwt.m.
%
%   Inputs:
%       order: Order of the fitting polynomials.
%       a: Start of the integration interval.
%       b: End of the integration interval.
%
%   Outputs:
%       points: Locations at which the integrand is evaluated.
%       weights: Corresponding weights for the quadrature sum.

    [points, weights] = external.lgwt.lgwt(order, a, b);
end
