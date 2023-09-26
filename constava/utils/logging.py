import logging
import logging.config

logger_configuration = {
    "version" : 1,
    "formatters" : {
        "default": {
            "format" : "[{asctime}] {message}",
            "datefmt" : "%Y-%m-%d %H:%M:%S",
            "style" : "{",
            "validate" : True
        },
    },
    "handlers" : {
        "console" : {
            "class" : "logging.StreamHandler",
            "formatter" : "default",
            "stream": "ext://sys.stdout",
        },
        "null" : {
            "class" : "logging.NullHandler",
        }
    },
    "loggers" : {
        "Constava" : {
            "handlers" : ["console"],
            "level" : "WARNING",
        },
        "Dummy" : {
            "handlers" : ["null"],
        },
    },
}

logging.config.dictConfig(logger_configuration)
