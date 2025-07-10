function mesh = filterMeshByGenome(originalMesh, genome)
% filterMeshByGenome: Returns a new mesh containing only the faces of the original mesh indicated by the genome mask.
%
%   Inputs:
%     originalMesh: Mesh object containing vertices and faces
%     genome: logical 1xF array, where F is the number of faces in originalMesh
%
%   Output:
%     mesh: New Mesh object with the same vertices and only the faces selected by the genome mask

    faceCount = size(originalMesh.faces(), 2);

    % Validate genome length
    if numel(genome) ~= faceCount
        error('Genome length must match the number of faces in originalMesh.');
    end

    % Select faces according to the genome mask
    allFaces = originalMesh.faces();
    selectedFaces = allFaces(:, genome);

    % Create a new mesh with selected faces
    mesh = rooftopMoM.mesh.Mesh(originalMesh.vertices(), selectedFaces);
end
