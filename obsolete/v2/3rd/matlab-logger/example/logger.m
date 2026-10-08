% Configure logger: console + file, debug level.
logger.Logger.configure( ...
    'minLevel', logger.LogLevel.DEBUG, ...
    'sinks', { ...
        logger.sink.ConsoleSink(), ...
        logger.sink.FileSink('app.log') ...
    } ...
);

% Use full API.
logger.Logger.info('Starting app: pid=%d', feature('getpid'));
logger.Logger.debug('Matrix sample: %s', mat2str(rand(1,3)));

% Use short-form wrappers.
INFO('Startup complete');
DEBUG('Debugging value x=%f', 3.14);
WARN('This is a warning - config missing, using defaults');
ERROR('This is an error (but not throwing)');

% Turn off logging quickly.
logger.Logger.setEnabled(false);
INFO('This will not be printed');

% Reconfigure the logger.
logger.Logger.configure( ...
    'timeFormat', 'ss:SSS', ...
    'template', 'time: %s level: %s message: %s\n' ...
);

% Turn on again.
logger.Logger.setEnabled(true);

INFO('Test message in the new format.');
