classdef QuadPoints
    properties
        order_ uint32;
        points2Du_ double;
        points2Dv_ double;
        weights4D_ double;
    end

    methods
        function obj = QuadPoints(order)
            function [points, weights] = generateQuadPoints1D(order)
                [points, weights] = mom.lgwt(double(order), 0, 1);
            end

            [p1D, w1D] = generateQuadPoints1D(order);

            [p2Du, p2Dv] = ndgrid(p1D, p1D);
            obj.points2Du_ = p2Du(:).';
            obj.points2Dv_ = p2Dv(:).';

            [w4D1, w4D2, w4D3, w4D4] = ndgrid(w1D, w1D, w1D, w1D);
            obj.weights4D_ = w4D1(:) .* w4D2(:) .* w4D3(:) .* w4D4(:);
        end
    end

    methods
        function points = getTransformed(obj, v1, v2, v4)
            uSide= v2 - v1;
            vSide = v4 - v1;

            points = uSide .* obj.points2Du_ + vSide .* obj.points2Dv_ + v1;
        end
    end
end