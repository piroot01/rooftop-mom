clc; clearvars;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', '..', 'src');
addpath(srcDir);

loggerDir = fullfile(currentDir, '..', '..', '3rd', 'matlab-logger', 'src');
addpath(loggerDir);

% Configure the logger.
logger.Logger.configure( ...
    'timeFormat', 'ss.SSS', ...
    'minLevel', logger.LogLevel.TRACE, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
    } ...
);

countX = 10;
countY = 10;
sizeX = 1;
sizeY = 1;

order = 3;
elementCount = 1e6;
k = pi / (sizeY * countY);
V = 1;

quad = rooftopMoM.basisFunctions.GaussLegendre2RectQuadrature(order);
rectPrimitive = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);

rooftops = rooftopMoM.basisFunctions.UniformRectPrimitiveUniformRooftops(rectPrimitive, quad);
rooftops.generate();

planner = rooftopMoM.core.greensMatrix.AssemblerChunkPlanner(elementCount);
context = rooftopMoM.core.impedanceMatrix.AssemblerContext(quad, rooftops, k);
assembler = rooftopMoM.core.greensMatrix.FullAssembler(planner);

greensMatrix = assembler.generate(context);

spectralGreensMatrix = greensMatrix.generateSpectralGreensMatrix();

greensMatrixApplicator = rooftopMoM.core.greensMatrix.FFTAcceleratedApplicator( ...
    context, spectralGreensMatrix, greensMatrix.rooftopContent());

preconditionerApplicator = rooftopMoM.core.greensMatrix.preconditioner.FFTAcceleratedApplicator( ...
    context, spectralGreensMatrix, greensMatrix.rooftopContent());
preconditionerApplicator.generate();

solver = rooftopMoM.core.solver.PreconditionedConjugateGradient(greensMatrixApplicator, preconditionerApplicator);

mode = rooftopMoM.operator.utils.generateMiddlePortMode(rooftops.count());
excitation = rooftopMoM.operator.DeltaGapFeederExcitation.generate(rooftops.count(), mode, V);

tolerance = 1e-3;
maxIterCount = 1e3;

outputConvSolver = solver.solve( ...
    'excitation', excitation, ...
    'tolerance', tolerance, ...
    'maxIterCount', maxIterCount ...
);

fprintf('Number of the iterations: %d\n', outputConvSolver.count());

currentSpectral = outputConvSolver.postRunInfo().current();

%% Caledron preconditioner setup
[circSizeX, circSizeY] = size(spectralGreensMatrix.blockXX());

[u, v] = ndgrid(0:(circSizeX - 1), 0:(circSizeY - 1));

phaseHalfX = exp((-1i * pi * u) / circSizeX);
phaseHalfY = exp((-1i * pi * v) / circSizeY);

% Gram matrix
gramStencil = [1 6 1; 6 36 6; 1 6 1] / 64;
gramPad = zeros(circSizeX, circSizeY);
gramPad(1:3, 1:3) = gramStencil;
gramPad = circshift(gramPad, [-1 -1]);
gramSpectr = fft2(gramPad);

invGramSpectr = 1 ./ (gramSpectr + 1e-12);

applyHandle = @(x) applyCalderonOperator(x, greensMatrixApplicator, phaseHalfX, phaseHalfY, invGramSpectr, rooftops, spectralGreensMatrix, true);

rhs = applyCalderonOperator(excitation, greensMatrixApplicator, phaseHalfX, phaseHalfY, invGramSpectr, rooftops, spectralGreensMatrix, false);

[currentCaledron, flag, ~, iterCount] = cgs(applyHandle, rhs, tolerance, maxIterCount);

norm(currentCaledron - currentSpectral)

iterCount

function output = applyCalderonOperator(input, greensMatrixApplicator, phaseHalfX, phaseHalfY, invGramSpectr, rooftops, spectralMatrixData, applyPrimal)
    % wrapper for the full A^2 preconditioned operator
    
    % =================================================================
    % PASS 1: Primal Grid (J_primal -> E_primal)
    % =================================================================
    % This calls your existing code exactly as written
    if applyPrimal
        E_primal = greensMatrixApplicator.apply(input);
    else
        E_primal = input;
    end
    
    
    % =================================================================
    % THE LINK: Project Primal Field -> Dual Current
    % Operation: J_dual = IFFT( FFT(E_primal) * Shift * G^-1 )
    % =================================================================
    
    % 1. Recover structure from the vector (Reshape)
    %    (We need to separate X and Y components to apply shifts correctly)
    primitive = rooftops.primitive();
    basisSizeX = rooftops.countDir('x');
    
    % Reshape X-component
    Ex_reshaped = reshape(E_primal(1:basisSizeX), primitive.countX() - 1, primitive.countY());
    % Reshape Y-component
    Ey_reshaped = reshape(E_primal((basisSizeX + 1):end), primitive.countX(), primitive.countY() - 1);
    
    % 2. Move to Spectral Domain
    [blockSizeX, blockSizeY] = spectralMatrixData.blockSize();
    Ex_spectr = fft2(Ex_reshaped, blockSizeX, blockSizeY);
    Ey_spectr = fft2(Ey_reshaped, blockSizeX, blockSizeY);
    
    % 3. Apply SHIFT and INVERSE GRAM
    %    We shift the "Voltage" testing location to the "Current" source location on the dual grid.
    Jx_dual_spectr = Ex_spectr .* phaseHalfX .* phaseHalfY .* invGramSpectr;
    Jy_dual_spectr = Ey_spectr .* phaseHalfX .* phaseHalfY .* invGramSpectr;
    
    % 4. Back to Spatial Domain (Truncate to metal shape)
    %    Crucial: This truncation forces the dual current to respect the metal boundaries.
    Jx_dual = ifft2(Jx_dual_spectr);
    Jy_dual = ifft2(Jy_dual_spectr);
    
    % Truncate (Dual grid assumes same number of unknowns for simplicity in this mapping)
    Jx_dual_trunc = Jx_dual(1:(primitive.countX() - 1), 1:primitive.countY());
    Jy_dual_trunc = Jy_dual(1:primitive.countX(), 1:(primitive.countY() - 1));
    
    % Re-pack into vector for the next pass
    J_dual_vector = [Jx_dual_trunc(:); Jy_dual_trunc(:)];
    
    
    % =================================================================
    % PASS 2: Dual Grid (J_dual -> E_dual)
    % =================================================================
    % We reuse the EXACT same Z-matrix logic.
    % The physics is translation invariant, so Z is the same.
    E_dual = greensMatrixApplicator.apply(J_dual_vector);
    
    
    % =================================================================
    % UN-SHIFT: Align back to Primal Grid
    % Operation: V_out = IFFT( FFT(E_dual) * Shift_Conj )
    % =================================================================
    
    % Unpack
    E_dual_x = reshape(E_dual(1:basisSizeX), primitive.countX() - 1, primitive.countY());
    E_dual_y = reshape(E_dual((basisSizeX + 1):end), primitive.countX(), primitive.countY() - 1);
    
    % Spectral
    Ed_x_spectr = fft2(E_dual_x, blockSizeX, blockSizeY);
    Ed_y_spectr = fft2(E_dual_y, blockSizeX, blockSizeY);
    
    % Apply Conjugate Shift (Move back by -Delta/2)
    V_out_x_spectr = Ed_x_spectr .* conj(phaseHalfX) .* conj(phaseHalfY);
    V_out_y_spectr = Ed_y_spectr .* conj(phaseHalfX) .* conj(phaseHalfY);
    
    % IFFT and Final Truncation
    V_out_x = ifft2(V_out_x_spectr);
    V_out_y = ifft2(V_out_y_spectr);
    
    V_final_x = V_out_x(1:(primitive.countX() - 1), 1:primitive.countY());
    V_final_y = V_out_y(1:primitive.countX(), 1:(primitive.countY() - 1));
    
    % Final Output Vector
    output = [V_final_x(:); V_final_y(:)];
    
end

%{
                    primitive = obj.context_.rooftops().primitive();
                    countX = primitive.countX();
                    countY = primitive.countY();

                    % First we process the x direction of the input
                    inputXReshaped = reshape(input(1:basisSizeX), countX - 1, countY);

                    [blockSizeX, blockSizeY] = obj.spectralMatrixData_.blockSize();
                    inputXSpectr = fft2(inputXReshaped, blockSizeX, blockSizeY);
                    clear inputXReshaped;

                    outputXSpectr = obj.spectralMatrixData_.blockXX() .* inputXSpectr;
                    outputYSpectr = obj.spectralMatrixData_.blockYX() .* inputXSpectr;
                    clear inputXSpectr;

                    % Secondly we process the y direction of the input
                    inputYReshaped = reshape(input((basisSizeX + 1):end), countX, countY - 1);

                    inputYSpectr = fft2(inputYReshaped, blockSizeX, blockSizeY);
                    clear inputYReshaped;

                    outputXSpectr = outputXSpectr + obj.spectralMatrixData_.blockXY() .* inputYSpectr;
                    outputYSpectr = outputYSpectr + obj.spectralMatrixData_.blockYY() .* inputYSpectr;
                    clear inputYSpectr;

                    % And apply the inverse FFTs sequentially
                    outputX = ifft2(outputXSpectr);
                    clear outputXSpectr; % spatial domain

                    outputXTruncated = outputX(1:(countX - 1), 1:countY);
                    clear outputX;

                    outputY = ifft2(outputYSpectr);
                    clear outputYSpectr;

                    outputYTruncated = outputY(1:countX, 1:(countY - 1));
                    clear outputY;

                    % final concatenated output
                    output = [outputXTruncated(:); outputYTruncated(:)];
%}
