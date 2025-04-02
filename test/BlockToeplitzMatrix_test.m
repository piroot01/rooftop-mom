clc; clearvars; format long;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

countX = 3;
countY = 3;
sizeX = 1;
sizeY = 1;
origin = [0, 0, 0];
order = 2;

lambda = 5;
k = 2 * pi / lambda;

memoryElementwise = 16;
memory_MB = memoryElementwise * 8 / 1e6;

% 1. step: generate a mesh
myRectMeshPrimitive1 = mom.MeshRectPrimitive(countX, countY, sizeX, sizeY, []);
myRectMeshPrimitive2 = mom.MeshRectPrimitive(countX, countY, sizeX * 1.23, sizeY, []);
myRectMesh1 = mom.MeshRect(myRectMeshPrimitive1, origin, 0, 0, 0);
myRectMesh2 = mom.MeshRect(myRectMeshPrimitive2, origin + [0, 0, 1], 0, 0, 0);
myMeshCluster = mom.Mesh2DCluster();
myMeshCluster.loadMesh(myRectMesh1.getAsMesh());
%myMeshCluster.loadMesh(myRectMesh2.getAsMesh());

% 2. step: create Quadrature
quadrature = mom.Quadrature(order);

% 3. step: create BasisFunctions
basisFuncs = mom.BasisFunctions(myMeshCluster);

basisFuncs.getFaces()

% 4. step: create RooftopBasisFunctions
rooftops = mom.RooftopBasisFunctions(myMeshCluster, basisFuncs.getFaces());

% 5. step: generate quadrature points and rooftop points
rooftops.generatePoints(quadrature);

stratCtx = mom.ImpedanceMatrix.StrategyContext(quadrature, basisFuncs, rooftops, k);
compStratCtx = mom.ImpedanceMatrix.ComputeChunkedContext(memory_MB);
%combStrat = mom.ImpedanceMatrix.ComputeChunked(compStratCtx);
compStrat = mom.ImpedanceMatrix.ComputeMagPartChunked(compStratCtx);

Z = compStrat.compute(stratCtx);

disp(Z);

% magnetical part of the impedance matrix
% we will compute the interaction between first and third basis function
%allQuadraturePoints = rooftops.getQuadraturePoints()
allRooftopPoints = rooftops.getRooftopPoints();

quadratureWeights = quadrature.getWeights();

function points = getQuadraturePoints(quadrature, m)
    points = quadrature.getQuadraturePoints();
    points = points + [(m - [1, 1]).'; 0];
end

% obtain the index of the first basisFunction in the y direction
facesX = basisFuncs.getFacesXAsCell();
basisFuncsXCount = numel(facesX{:}) / 2;
basisFunctionsYStartingIndex = basisFuncsXCount;

function greensFunction = computeGreensFunction(k, quadraturePoints0, quadraturePointsShifted)
    distanceMatrix = mom.DistanceMatrix.getChunk(quadraturePoints0, quadraturePointsShifted);
    distanceMatrixZero = (distanceMatrix < eps);
    greensFunction = exp(-1i * k * distanceMatrix) ./ distanceMatrix;
    greensFunction(distanceMatrixZero) = 0;
end

function result = compute_hSemi(greensFunction, rooftopPointsAlpha, rooftopPointsBeta, weights)
    rooftopDot = rooftopPointsAlpha.' * rooftopPointsBeta;
    integrand = greensFunction .* rooftopDot;
    integral = integrand .* weights;
    result = sum(integral, 'all');
end

function h1 = compute_h1(k, m, n, alpha, beta, quadrature, rooftops, weights, gridSize)
    if (beta == 'x' && m(1) == gridSize(1)) || (beta == 'y' && m(2) == gridSize(2))
        h1 = 0;
    else
        quadraturePointsAlphaPlus = getQuadraturePoints(quadrature, n);
        quadraturePoints0 = getQuadraturePoints(quadrature, [1, 1]);
        rooftopPointsAlphaPlus = rooftops.getWedgePoints(quadraturePoints0, alpha, 'p');
        rooftopPointsAlphaMinus = rooftops.getWedgePoints(quadraturePoints0, alpha, 'm');

        %fprintf("[m, n, alpha, beta] = [%s, %s, %s, %s]\n", mat2str(m), mat2str(n), alpha, beta);
        quadraturePointsBetaPlus = getQuadraturePoints(quadrature, m);
        rooftopPointsBetaPlus = rooftops.getWedgePoints(quadraturePoints0, beta, 'p');

        greensFunctionPlus = computeGreensFunction(k, quadraturePointsAlphaPlus, quadraturePointsBetaPlus);

        if alpha == 'x'
            quadraturePointsAlphaPlus = quadraturePointsAlphaPlus + [1, 0, 0].';
        else
            quadraturePointsAlphaPlus = quadraturePointsAlphaPlus + [0, 1, 0].';
        end

        greensFunctionMinus = computeGreensFunction(k, quadraturePointsAlphaPlus, quadraturePointsBetaPlus);


        h1Plus = compute_hSemi(greensFunctionPlus, rooftopPointsAlphaPlus, rooftopPointsBetaPlus, weights);
        h1Minus = compute_hSemi(greensFunctionMinus, rooftopPointsAlphaMinus, rooftopPointsBetaPlus, weights);

        h1 = h1Plus + h1Minus;
    end
end

function h2 = compute_h2(k, m, n, alpha, beta, quadrature, rooftops, weights)
    if (beta == 'x' && m(1) == 1) || (beta == 'y' && m(2) == 1)
        h2 = 0;
    else
        %fprintf("[m, n, alpha, beta] = [%s, %s, %s, %s]\n", mat2str(m), mat2str(n), alpha, beta);
        quadraturePointsAlphaPlus = getQuadraturePoints(quadrature, n);
        quadraturePoints0 = getQuadraturePoints(quadrature, [1, 1]);
        rooftopPointsAlphaPlus = rooftops.getWedgePoints(quadraturePoints0, alpha, 'p');
        rooftopPointsAlphaMinus = rooftops.getWedgePoints(quadraturePoints0, alpha, 'm');

        quadraturePointsBetaMinus = getQuadraturePoints(quadrature, m);
        rooftopPointsBetaMinus = rooftops.getWedgePoints(quadraturePoints0, beta, 'm');

        greensFunctionPlus = computeGreensFunction(k, quadraturePointsAlphaPlus, quadraturePointsBetaMinus);

        if alpha == 'x'
            quadraturePointsAlphaPlus = quadraturePointsAlphaPlus + [1, 0, 0].';
        else
            quadraturePointsAlphaPlus = quadraturePointsAlphaPlus + [0, 1, 0].';
        end

        greensFunctionMinus = computeGreensFunction(k, quadraturePointsAlphaPlus, quadraturePointsBetaMinus);

        h2Plus = compute_hSemi(greensFunctionPlus, rooftopPointsAlphaPlus, rooftopPointsBetaMinus, weights);
        h2Minus = compute_hSemi(greensFunctionMinus, rooftopPointsAlphaMinus, rooftopPointsBetaMinus, weights);

        h2 = h2Plus + h2Minus;
    end
end

%h1 = compute_h1(k, [1, 2], 'x', 'x', quadrature, rooftops, quadratureWeights, [2, 2]);
%h2 = compute_h2(k, [2, 2], 'x', 'x', quadrature, rooftops, quadratureWeights);

gridSize = [3, 3];

directions = ['x', 'y'];

G1 = zeros(3, 3, 2, 2);
G2 = zeros(3, 3, 2, 2);

n = [1, 1];
for my = 1:3
    for mx = 1:3
        m = [mx, my];
        H1 = zeros(2);
        H2 = zeros(2);
        for a = 1:2
            for b = 1:2
                alpha = directions(a);
                beta = directions(b);
                H1(a, b) = compute_h1(k, m, n, alpha, beta, quadrature, rooftops, quadratureWeights, gridSize);
                H2(a, b) = compute_h2(k, m, n, alpha, beta, quadrature, rooftops, quadratureWeights);
            end
        end

        G1(mx, my, :, :) = H1;
        G2(mx, my, :, :) = H2;
    end
end

%Z0 = 4 * pi * 29.9792458;

%G1 = (1 / (4 * pi)) * 1i * k * Z0 * G1
%G2 = (1 / (4 * pi)) * 1i * k * Z0 * G2

%G1(1, 2, :, :)
G2(2, 2, :, :)

n = [2, 1];
H1 = zeros(2);
for a = 1:2
    for b = 1:2
        alpha = directions(a);
        beta = directions(b);
        H1(a, b) = compute_h1(k, [3, 2], n, alpha, beta, quadrature, rooftops, quadratureWeights, gridSize);
        H2(a, b) = compute_h2(k, [3, 2], n, alpha, beta, quadrature, rooftops, quadratureWeights);
    end
end

%H1
H2

%G1(1, 1, :, 1) + G2(2, 1, :, 1)
%G1(1, 1, :, 2) + G2(1, 2, :, 2)
%G1(1, 2, :, 1) + G2(2, 2, :, 1)
%G1(2, 1, :, 2) + G2(2, 2, :, 2)

%disp((1 / (4 * pi)) * 1i * k * Z0 * (h1 + h2));
