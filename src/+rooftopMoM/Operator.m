classdef Operator < rooftopMoM.Product
% Operator: Matrix-free product with the Green's matrix of a Problem, or with a part of it.
%       operator = rooftopMoM.Operator(problem, 'product', 'splitFFT');
%       output = operator.apply(input);

    methods
        function obj = Operator(problem, options)
        % Operator: Creates the product; the problem assembles its generators if it has not yet.
        %
        %   Input:
        %       problem: rooftopMoM.Problem.
        %
        %   Optional name-value inputs:
        %       'product': How the product is computed (default 'splitFFT'):
        %           'fft': FFTs of the circulant embedding, about twice the grid size in every direction.
        %           'splitFFT': FFTs of the grid size only, the fastest from about 128 x 128 faces up.
        %           'leanSplitFFT': The split FFT in less working memory, but slower.
        %       'modifier': Function handle applied to every generator, e.g. @real or @imag for the real or the
        %           imaginary part of the Green's matrix (default none).

            arguments
                problem (1, 1) rooftopMoM.Problem
                options.product (1, :) char {mustBeMember(options.product, {'fft', 'splitFFT', 'leanSplitFFT'})} = ...
                    'splitFFT'
                options.modifier (1, 1) function_handle = @(generator) generator
            end

            generators = problem.generators().modify(options.modifier);
            rooftops = problem.rooftops();
            [sizeX, sizeY] = problem.mesh().faceSize();

            switch options.product
                case 'fft'
                    spectra = rooftopMoM.assembler.greensMatrix.utils.generateEmbeddingSpectra(generators, sizeX, ...
                        sizeY);
                    applicator = rooftopMoM.assembler.greensMatrix.FFTApplicator(rooftops, spectra, ...
                        generators.content());
                case 'splitFFT'
                    spectra = rooftopMoM.assembler.greensMatrix.utils.generateSubbandSpectra(generators, sizeX, sizeY);
                    applicator = rooftopMoM.assembler.greensMatrix.SplitFFTApplicator(rooftops, spectra, ...
                        generators.content());
                case 'leanSplitFFT'
                    spectra = rooftopMoM.assembler.greensMatrix.utils.generateSubbandSpectra(generators, sizeX, sizeY);
                    applicator = rooftopMoM.assembler.greensMatrix.LeanSplitFFTApplicator(rooftops, spectra, ...
                        generators.content());
            end

            obj.problem_ = problem;
            obj.applicator_ = applicator;
        end
    end
end
