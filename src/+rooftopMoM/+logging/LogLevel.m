classdef LogLevel < uint8
% LogLevel: Levels of the log messages, from the least to the most verbose.
%   The logger writes the messages up to its own level. OFF writes nothing.
%       ERROR: failures that stop the computation.
%       WARNING: results to check, e.g. a solver that did not converge (default level).
%       INFO: one line per stage, with its time and memory.
%       DEBUG: inputs, chunk progress and convergence.
%       TRACE: one line per iteration and per chunk.

    enumeration
        OFF         (0)
        ERROR       (1)
        WARNING     (2)
        INFO        (3)
        DEBUG       (4)
        TRACE       (5)
    end
end
