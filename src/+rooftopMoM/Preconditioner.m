classdef Preconditioner < rooftopMoM.Product
% Preconditioner: Circulant preconditioner of the Green's matrix of a Problem, the product with the inverse of a block
%   circulant approximation of the matrix.
%       preconditioner = rooftopMoM.Preconditioner(problem, 'kernel', 'tChan', 'builder', 'partwise');

    methods
        function obj = Preconditioner(problem, options)
        % Preconditioner: Creates the preconditioner; the problem assembles its generators if it has not yet.
        %
        %   Input:
        %       problem: rooftopMoM.Problem.
        %
        %   Optional name-value inputs:
        %       'kernel': Circulant approximation of a Toeplitz generator, 'strang', 'tChan', 'rChan' or 'huckle'
        %           (default 'tChan').
        %       'builder': What is approximated (default 'partwise'):
        %           'partwise': Each of the three generators.
        %           'blockwise': Each of the four blocks of the matrix.
        %           'hybrid': The diagonal blocks and the electric generator.
        %       'bandwidth': 1x2 bandwidths of the 'huckle' kernel in faces along x and y (default the whole grid).
        %       'modifier': Function handle applied to every generator first, e.g. @imag for the imaginary part of
        %           the Green's matrix (default none).

            arguments
                problem (1, 1) rooftopMoM.Problem
                options.kernel (1, :) char {mustBeMember(options.kernel, {'strang', 'tChan', 'rChan', 'huckle'})} = ...
                    'tChan'
                options.builder (1, :) char {mustBeMember(options.builder, {'partwise', 'blockwise', 'hybrid'})} = ...
                    'partwise'
                options.bandwidth {mustBeInteger, mustBePositive, mustBePairOrEmpty_} = []
                options.modifier (1, 1) function_handle = @(generator) generator
            end

            bandwidth = options.bandwidth;

            if isempty(bandwidth)
                [countX, countY] = problem.mesh().count();
                bandwidth = [countX, countY];
            end

            kernels = struct( ...
                'strang', @rooftopMoM.assembler.preconditioner.kernels.generateStrang, ...
                'tChan', @rooftopMoM.assembler.preconditioner.kernels.generateTChan, ...
                'rChan', @rooftopMoM.assembler.preconditioner.kernels.generateRChan, ...
                'huckle', @(generator) rooftopMoM.assembler.preconditioner.kernels.generateHuckle(generator, ...
                    bandwidth(1), bandwidth(2)));
            builders = struct( ...
                'partwise', @rooftopMoM.assembler.preconditioner.generatePartwiseSpectra, ...
                'blockwise', @rooftopMoM.assembler.preconditioner.generateBlockwiseSpectra, ...
                'hybrid', @rooftopMoM.assembler.preconditioner.generateHybridSpectra);

            obj.problem_ = problem;
            obj.applicator_ = rooftopMoM.assembler.preconditioner.generateCirculant(problem.rooftops(), ...
                problem.generators().modify(options.modifier), kernels.(options.kernel), builders.(options.builder));

            % A zero determinant at one frequency leaves Inf or NaN in the inverse.
            spectra = obj.applicator_.spectra();
            blocks = {spectra.blockXX(), spectra.blockXY(), spectra.blockYX(), spectra.blockYY()};

            if ~all(cellfun(@allfinite, blocks))
                error('rooftopMoM:preconditioner:singular', ['The circulant approximation is singular at some ' ...
                    'frequency. Use another kernel or builder.']);
            end
        end
    end
end

function mustBePairOrEmpty_(value)
    if ~isempty(value) && numel(value) ~= 2
        error('rooftopMoM:preconditioner:invalidBandwidth', 'The bandwidth is a pair [bandwidthX, bandwidthY].');
    end
end
