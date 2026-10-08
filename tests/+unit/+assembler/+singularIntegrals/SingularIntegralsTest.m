classdef SingularIntegralsTest < unit.Base
% SingularIntegralsTest: Tests of rooftopMoM.assembler.singularIntegrals against adaptive integration.
%   The four-dimensional integrals depend only on the differences nu = u - u' and mu = v - v'. Over the unit square,
%   a weight w(u, u') becomes the correlation W(nu), and the integrals reduce to
%       integral over [0, 1]^2 of (W(nu) + W(-nu)) * 2 * (1 - mu) / R(nu, mu),
%   which integral2 evaluates accurately. The correlations are
%       1: W(nu) + W(-nu) = 2 * (1 - nu);
%       u * u': W(nu) + W(-nu) = (1 - nu)^2 * (2 + nu) / 3;
%       u * (1 - u'): W(nu) + W(-nu) = (1 - nu^3) / 3.

    properties (TestParameter)
        sizes = {[0.1, 0.15], [1, 1], [0.02, 0.1]};
    end

    methods (Test)
        function testType1(testCase, sizes)
            expected = testCase.reducedIntegral_(sizes, @(nu) (1 - nu).^2 .* (2 + nu) / 3);
            testCase.verifyEqual(rooftopMoM.assembler.singularIntegrals.type1(sizes(1), sizes(2)), expected, ...
                'RelTol', 1e-10);
        end

        function testType2(testCase, sizes)
            expected = testCase.reducedIntegral_(sizes, @(nu) (1 - nu.^3) / 3);
            testCase.verifyEqual(rooftopMoM.assembler.singularIntegrals.type2(sizes(1), sizes(2)), expected, ...
                'RelTol', 1e-10);
        end

        function testType3(testCase, sizes)
            expected = testCase.reducedIntegral_(sizes, @(nu) 2 * (1 - nu));
            testCase.verifyEqual(rooftopMoM.assembler.singularIntegrals.type3(sizes(1), sizes(2)), expected, ...
                'RelTol', 1e-10);
        end
    end

    methods (Access = private)
        function value = reducedIntegral_(~, sizes, correlation)
            integrand = @(nu, mu) correlation(nu) .* 2 .* (1 - mu) ./ sqrt((sizes(1) * nu).^2 + (sizes(2) * mu).^2);
            value = integral2(integrand, 0, 1, 0, 1, 'Method', 'iterated', 'AbsTol', 1e-14, 'RelTol', 1e-12);
        end
    end
end
