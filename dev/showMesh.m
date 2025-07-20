clc; clearvars; close all;

currentDir = fileparts(mfilename('fullpath'));
srcDir = fullfile(currentDir, '..', 'src');
addpath(srcDir);

% First mesh
countX = 3;
countY = 3;
sizeX = 1;
sizeY = 1;

rectPrimitive1 = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect1 = rooftopMoM.mesh.UniformRect(rectPrimitive1);
mesh1 = rect1.generateMesh();

genome = true(countX * countY, 1);
genome([2, 4, 8, 6]) = false;

localVertices = projectVerticesToLocalCoords(mesh1.vertices{1}, mesh1.faces{1});

displayMesh(mesh1.faces{1}, localVertices, genome);

% Second mesh
countX = 20;
countY = 20;
sizeX = 0.15;
sizeY = 0.15;

rectPrimitive2 = rooftopMoM.mesh.UniformRectPrimitive(countX, countY, sizeX, sizeY);
rect2 = rooftopMoM.mesh.UniformRect(rectPrimitive2);
mesh2 = rect2.generateMesh();

load('genome.mat');

threshold = 0.01;

genome = logical(newGenome.' > threshold);

localVertices = projectVerticesToLocalCoords(mesh2.vertices{1}, mesh2.faces{1});

displayMesh(mesh2.faces{1}, localVertices, genome);

function points = projectVerticesToLocalCoords(vertices, faces, faceIndex)
    if nargin < 3 || isempty(faceIndex)
        faceIndex = 1;
    end

    % Compute local basis vectors.
    [u, v] = rooftopMoM.mesh.utils.generateRectLocalCoords(vertices, faces, faceIndex);

    points = zeros(2, size(vertices, 2));
    % Project vertices onto local axes.
    points(1, :) = u.' * vertices;
    points(2, :) = v.' * vertices;
end

function displayMesh(faces, verts, genome)
%DISPLAYMESH Display a 2D quad mesh and color faces by a binary genome
%
%   displayMesh(faces, verts, genome)
%       faces : 4xN matrix, each column holds indices of the 4 vertices of a quad
%       verts : 2xM matrix, each column is the (x; y) coordinates of a vertex
%       genome: 1xN or Nx1 logical/binary array, face active (1) -> red, inactive (0) -> white

    % Input validation
    [fRows, numFaces] = size(faces);
    [vRows, ~] = size(verts);
    if fRows ~= 4
        error('faces must be a 4xN matrix');
    end
    if vRows ~= 2
        error('verts must be a 2xM matrix');
    end
    if numel(genome) ~= numFaces
        error('genome length must match number of faces');
    end
    genome = logical(genome(:));  % ensure logical column vector

    figure;
    hold on;

    % Loop through each face and draw it
    for i = 1:numFaces
        vidx = faces(:, i);
        xy = verts(:, vidx);  % 2x4
        % Choose color based on genome
        if genome(i)
            fc = [1, 0, 0];  % red for active
        else
            fc = [1, 1, 1];  % white for inactive
        end
        patch('XData', xy(1, :), 'YData', xy(2, :), ...
              'FaceColor', fc, 'EdgeColor', [0,0,0]);
    end

    axis equal;
    axis off;
    hold off;
end