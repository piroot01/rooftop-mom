classdef AssemblerContext < handle
    properties (Access = private)
        quadrature_;
        rooftops_;

        % 1xM vector of wavenumbers in which the Assember will evaluate the impedance matrix.
        wavenumber_;
    end

    methods
        function obj = AssemblerContext(quadrature, rooftops, wavenumber)
            obj.quadrature_ = quadrature;
            obj.rooftops_ = rooftops;
            obj.wavenumber_ = wavenumber;
        end
    end

    methods
        function quadrature = quadrature(obj)
            quadrature = obj.quadrature_;
        end

        function rooftops = rooftops(obj)
            rooftops = obj.rooftops_;
        end

        function wavenumber = wavenumber(obj)
            wavenumber = obj.wavenumber_;
        end
    end
end
