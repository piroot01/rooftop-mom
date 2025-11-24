classdef PlaneWaveExcitation < rooftopMoM.operator.OperatorInterface
    methods (Static)
        function V = generateWithAssemblerContext(context, direction, polarization, intensity, useCachedData)
            if nargin < 5
                useCachedData = false;
            end

            V = rooftopMoM.operator.PlaneWaveExcitation.generate( ...
                context.quadrature(), ...
                context.rooftops(), ...
                context.wavenumber(), ...
                direction, ...
                polarization, ...
                intensity, ...
                useCachedData ...
            );
        end

        function V = generate(quadrature, rooftops, wavenumber, direction, polarization, intensity, useCachedData)
            if nargin < 7
                useCachedData = false;
            end

            [quadraturePoints, rooftopPoints, rooftopValues] = rooftopMoM.operator.utils.generateRooftopData(rooftops, useCachedData);

            direction = direction(:).';
            polarization = polarization(:).';

            matrixCount = numel(wavenumber);
            basisSize = rooftops.count();
            order = quadrature.order();
            blockSize = order^2;

            weights = quadrature.weights2D();

            rooftopEdgeLengths = rooftopMoM.operator.utils.generateRooftopEdgeLengths(rooftops, rooftopValues);

            computeContext = struct( ...
                'basisSize', basisSize, ...
                'blockSize', blockSize, ...
                'quadraturePoints', quadraturePoints, ...
                'rooftopPoints', rooftopPoints, ...
                'rooftopEdgeLengths', rooftopEdgeLengths, ...
                'weights', repmat(weights, basisSize, 1).'  ...
            );

            V = complex(zeros(basisSize, 1, matrixCount));

            for n = 1:matrixCount
                V(:, :, n) = rooftopMoM.operator.PlaneWaveExcitation.generateImpl_(computeContext, direction, polarization, intensity, wavenumber(n));
            end
        end
    end

    methods (Static, Access = private)
        function V = generateImpl_(computeContext, direction, polarization, intensity, wavenumber)
            exponent = 1i * wavenumber * pagemtimes(direction, computeContext.quadraturePoints);

            integrand = computeContext.rooftopPoints .* exp(exponent);

            integral = pagemtimes(polarization, (integrand .* computeContext.weights));

            integralReshaped = reshape(integral, computeContext.blockSize, computeContext.basisSize, 2);

            V = sum(sum(integralReshaped, 1), 3) .* computeContext.rooftopEdgeLengths;

            V = intensity * V;
        end
    end
end
