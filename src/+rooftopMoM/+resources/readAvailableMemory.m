function bytes = readAvailableMemory()
% readAvailableMemory: Returns the memory the machine can give to new data without swapping. Linux reports it in
%   /proc/meminfo, Windows through the Java runtime of MATLAB (without it through the function memory) and macOS
%   through vm_stat.
%
%   Output:
%       bytes: Available memory in bytes, NaN if it cannot be read.

    bytes = NaN;

    try
        if ispc && usejava('jvm')
            % The same value as the function memory gives, in about 0.1 ms instead of 46 ms.
            bytes = java.lang.management.ManagementFactory.getOperatingSystemMXBean().getFreePhysicalMemorySize();
        elseif ispc
            [~, systemView] = memory();
            bytes = systemView.PhysicalMemory.Available;
        elseif ismac
            [status, text] = system('vm_stat');

            if status == 0
                bytes = rooftopMoM.resources.utils.parseVmStat(text);
            end
        else
            bytes = rooftopMoM.resources.utils.parseMeminfo(fileread('/proc/meminfo'));
        end
    catch exception
        rooftopMoM.logging.Logger.instance().warning('Cannot read the available memory: %s', exception.message);
    end
end
