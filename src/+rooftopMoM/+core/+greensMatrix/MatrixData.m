classdef MatrixData < handle
    properties (Access = private)
        data_ struct;
    end

    methods
        function obj = MatrixData(varargin)
            p = inputParser();

            addParameter(p, 'gxxFFT', []);
            addParameter(p, 'gxyFFT', []);
            addParameter(p, 'gyxFFT', []);
            addParameter(p, 'gyyFFT', []);

            parse(p, varargin{:});

            obj.data_ = struct( ...
                'gxxFFT', p.Results.gxxFFT, ...
                'gxyFFT', p.Results.gxyFFT, ...
                'gyxFFT', p.Results.gyxFFT, ...
                'gyyFFT', p.Results.gyyFFT ...
            );
        end
    end

    methods
        function data = data(obj)
            data = obj.data_;
        end

        function gxxFFT = gxxFFT(obj)
            gxxFFT = obj.data_.gxxFFT;
        end

        function gxyFFT = gxyFFT(obj)
            gxyFFT = obj.data_.gxyFFT;
        end

        function gyxFFT = gyxFFT(obj)
            gyxFFT = obj.data_.gyxFFT;
        end

        function gyyFFT = gyyFFT(obj)
            gyyFFT = obj.data_.gyyFFT;
        end
    end
end
