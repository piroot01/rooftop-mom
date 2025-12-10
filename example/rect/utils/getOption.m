function val = getOption(opts, name, default)
    if isfield(opts, name)
        val = opts.(name);
    else
        val = default;
    end
end
