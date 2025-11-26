clc; clearvars;

countX = 2000;
countY = 2000;
order = 2;
repeatCount = 5;

count = (countX - 1) * countY + countX * (countY - 1);
orderSquared = order^2;

singleQuadraturePointsPositive = randn([3, orderSquared]);
singleQuadraturePointsNegative = randn([3, orderSquared]);
singleRooftopPointsPositive = randn([3, orderSquared]);
singleRooftopPointsNegative = randn([3, orderSquared]);

shifts = randn([3, count]);

cols = orderSquared * count;

for n = 1:repeatCount
    tic;
    %{
    quadraturePointsPositive = repmat(singleQuadraturePointsPositive, 1, count);
    quadraturePointsNegative = repmat(singleQuadraturePointsNegative, 1, count);
    shiftsDistributed = repelem(shifts, 1, orderSquared);
    rooftopPointsPositive = repmat(singleRooftopPointsPositive, 1, count) + shiftsDistributed;
    rooftopPointsNegative = repmat(singleRooftopPointsNegative, 1, count) + shiftsDistributed;

    quadraturePoints = cat(3, quadraturePointsPositive, quadraturePointsNegative);
    rooftopPoints = cat(3, rooftopPointsPositive, rooftopPointsNegative);
    %}
    
    quadraturePoints = zeros(3, cols, 2, 'like', singleRooftopPointsPositive);
    quadraturePoints(:, :, 1) = repmat(singleQuadraturePointsPositive, 1, count);
    quadraturePoints(:, :, 2) = repmat(singleQuadraturePointsNegative, 1, count);

    shiftsDistributed = repelem(shifts, 1, orderSquared);
    rooftopPoints = zeros(3, cols, 2, 'like', singleRooftopPointsPositive);
    rooftopPoints(:, :, 1) = repmat(singleRooftopPointsPositive, 1, count) + shiftsDistributed;
    rooftopPoints(:, :, 2) = repmat(singleRooftopPointsNegative, 1, count) + shiftsDistributed;
    
    %{
    quadraturePoints = zeros(3, cols, 2, 'like', singleRooftopPointsPositive);
    quadraturePoints(:, :, 1) = repmat(singleQuadraturePointsPositive, 1, count);
    quadraturePoints(:, :, 2) = repmat(singleQuadraturePointsNegative, 1, count);

    shiftsDistributed = repelem(shifts, 1, orderSquared);
    rooftopPoints = zeros(3, cols, 2, 'like', singleRooftopPointsPositive);
    rooftopPoints(:, :, 1) = repmat(singleRooftopPointsPositive, 1, count);
    rooftopPoints(:, :, 2) = repmat(singleRooftopPointsNegative, 1, count);
    rooftopPoints = rooftopPoints + shiftsDistributed;
    %}
    %{
    quadraturePoints = cat(3, repmat(singleQuadraturePointsPositive, 1, count), repmat(singleQuadraturePointsNegative, 1, count));

    shiftsDistributed = repelem(shifts, 1, orderSquared);
    rooftopPoints = cat(3, repmat(singleRooftopPointsPositive, 1, count), repmat(singleRooftopPointsNegative, 1, count)) + shiftsDistributed;
    %}
    toc;
end