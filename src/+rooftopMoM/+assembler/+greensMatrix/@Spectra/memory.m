function bytes = memory(obj)
% memory: Returns the bytes of the blocks of all elements of the (array) object.

    bytes = 0;

    for n = 1:numel(obj)
        bytes = bytes + rooftopMoM.assembler.utils.generateByteCount(obj(n).blockXX_, obj(n).blockXY_, ...
            obj(n).blockYX_, obj(n).blockYY_);
    end
end
