classdef BasisFunctions
    properties (Access = private)
        % store the rooftop basis functions as a pair of face indices
        functions_ uint32;
    end

    methods
        function obj = BasisFunctions(mesh)
            function edges = generateEdges(faces)
                faceCount = size(faces, 2);

                % each face has 4 edges, each edge will be represented with a pair of vertex indices and one face index
                edges = zeros(4 * faceCount, 3, 'uint32');

                m = 1;
                for n = 1:faceCount
                    face = faces(:, n);

                    edge1 = [min(face(1), face(2)), max(face(1), face(2))];
                    edge2 = [min(face(2), face(3)), max(face(2), face(3))];
                    edge3 = [min(face(3), face(4)), max(face(3), face(4))];
                    edge4 = [min(face(4), face(1)), max(face(4), face(1))];

                    edges(m, :) = [edge1, n];
                    edges(m + 1, :) = [edge2, n];
                    edges(m + 2, :) = [edge3, n];
                    edges(m + 3, :) = [edge4, n];

                    m = m + 4;
                end
            end

            function faces = getDoubleFaces(edges)
                sortedEdges = sortrows(edges, [1, 2]);
                diffEdges = diff(sortedEdges(:, 1:2));

                doubleEdgesIndices = find(all(diffEdges == 0, 2));

                face1 = sortedEdges(doubleEdgesIndices, 3);
                face2 = sortedEdges(doubleEdgesIndices + 1, 3);

                faces = [min(face1, face2), max(face1, face2)];
            end

            obj.functions_ = getDoubleFaces(generateEdges(mesh.getFaces())).';
        end

        function func = get(obj)
            func = obj.functions_;
        end
    end
end