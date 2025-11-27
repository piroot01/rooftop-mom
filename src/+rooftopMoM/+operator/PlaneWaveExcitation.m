classdef PlaneWaveExcitation < rooftopMoM.operator.OperatorInterface
% Currently works only for UniformRectPrimitiveUniformRooftops :C
    methods (Static)
        function V = generateWithAssemblerContext(context, direction, polarization, intensity, chunkSize)
            V = rooftopMoM.operator.PlaneWaveExcitation.generate( ...
                context.quadrature(), ...
                context.rooftops(), ...
                context.wavenumber(), ...
                direction, ...
                polarization, ...
                intensity, ...
                chunkSize ...
            );
        end

        function V = generate(quadrature, rooftops, wavenumber, direction, polarization, intensity, chunkSize)
            direction = direction(:).';
            polarization = polarization(:).';

            context = struct( ...
                'chunkSize', chunkSize, ...
                'quadrature', quadrature, ...
                'rooftops', rooftops ...
            );

            matrixCount = numel(wavenumber);
            basisSize = rooftops.count();

            V = complex(zeros(basisSize, 1, matrixCount));

            for n = 1:matrixCount
                V(:, :, n) = rooftopMoM.operator.PlaneWaveExcitation.generateImpl_(context, direction, polarization, intensity, wavenumber(n));
            end
        end

        function V = generatePerpendicular(rooftops, polarization, intensity)
            nx = rooftops.countDir('x');
            ny = rooftops.countDir('y');

            V = zeros(nx + ny, 1);

            V(1:nx) = polarization(1);
            V((nx + 1):end) = polarization(2);

            V = V * intensity;
        end
    end

    methods (Static, Access = private)
        function V = generateImpl_(context, direction, polarization, intensity, wavenumber)
            Vx = rooftopMoM.operator.PlaneWaveExcitation.generateDir_(context, direction, polarization, intensity, wavenumber, 'x');
            Vy = rooftopMoM.operator.PlaneWaveExcitation.generateDir_(context, direction, polarization, intensity, wavenumber, 'y');

            V = [Vx; Vy];
        end

        function V = generateDir_(context, direction, polarization, intensity, wavenumber, dir)
            quadrature = context.quadrature;
            rooftops = context.rooftops;
            primitive = rooftops.primitive();

            rooftopCount = rooftops.countDir(dir);

            if rooftopCount == 0
                V = [];
                return
            end

            orderSquared = quadrature.order()^2;

            [chunks, chunkCount] = rooftopMoM.operator.utils.generateChunks(rooftopCount, orderSquared, context.chunkSize);

            [rooftopPoints, rooftopValues] = context.rooftops.generateRooftopsDir(dir);

            V = complex(zeros(rooftops.countDir(dir), 1));

            for n = 1:chunkCount
                DEBUG("PlaneWaveExcitation.generateDir_: Assembling chunk %d/%d.", n, chunkCount);

                chunk = chunks{n};
                chunkIndices = chunk.indicesStart:chunk.indicesEnd;
                chunkSize = numel(chunkIndices);

                rooftopPointsPositive = repmat(rooftopPoints.rooftopPointsPositive, 1, chunkSize);
                rooftopPointsNegative = repmat(rooftopPoints.rooftopPointsNegative, 1, chunkSize);
                rooftopPointsCombined = cat(3, rooftopPointsPositive, rooftopPointsNegative);

                shifts = rooftopMoM.basisFunctions.utils.generateRectGridShiftsSliceDirection( ...
                    primitive.countX, primitive.sizeX, primitive.sizeY, ...
                    chunkIndices, ...
                    dir ...
                );

                quadraturePointsPositive = rooftopMoM.core.greensMatrix.utils.generateQuadraturePointsSlice(rooftopPoints.quadraturePointsPositive, shifts);
                quadraturePointsNegative = rooftopMoM.core.greensMatrix.utils.generateQuadraturePointsSlice(rooftopPoints.quadraturePointsNegative, shifts);
                quadraturePointsCombined = cat(3, quadraturePointsPositive, quadraturePointsNegative);

                context = struct( ...
                    'blockSize', orderSquared, ...
                    'chunkSize', chunkSize, ...
                    'quadraturePoints', quadraturePointsCombined, ...
                    'rooftopPoints', rooftopPointsCombined, ...
                    'rooftopEdgeLengths', repelem(rooftopValues.length, 1, chunkSize), ...
                    'weights', repmat(quadrature.weights2D().', 1, chunkSize) ...
                );

                V(chunkIndices) = rooftopMoM.operator.PlaneWaveExcitation.generateDirChunk_(context, direction, polarization, intensity, wavenumber);
            end
        end

        function V = generateDirChunk_(context, direction, polarization, intensity, wavenumber)
            exponent = 1i * wavenumber * pagemtimes(direction, context.quadraturePoints);

            integrand = context.rooftopPoints .* exp(exponent);

            integral = pagemtimes(polarization, (integrand .* context.weights));

            integralReshaped = reshape(integral, context.blockSize, context.chunkSize, 2);

            V = sum(sum(integralReshaped, 1), 3) .* context.rooftopEdgeLengths;

            V = intensity * V;
        end
    end
end
