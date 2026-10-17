"""Logging shared by the command-line tools."""

import logging
import sys

# Speed up logging for PyPy
logging._srcfile = None
logging.logThreads = False
logging.logProcesses = False
logging.logMultiprocessing = False

# Parent of the per-module loggers, so one set of handlers covers all of TACT
logger = logging.getLogger("tact")

_handlers: list[logging.Handler] = []


def setup_logging(output: str, verbose: int) -> None:
    """Send TACT's log messages to `{output}.log.txt` and to standard error.

    The log file gets INFO messages, or DEBUG messages with `-vv`. Standard error gets warnings,
    INFO messages with `-v`, or DEBUG messages with `-vv`. Calling this again replaces the handlers
    from the previous call.

    Args:
        output: Base name of the output files.
        verbose: How many times `-v` was given.
    """
    reset_logging()
    console_level = logging.DEBUG if verbose >= 2 else logging.INFO if verbose == 1 else logging.WARNING
    file_handler = logging.FileHandler(output + ".log.txt")
    file_handler.setLevel(min(console_level, logging.INFO))
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(console_level)
    for handler in (file_handler, console_handler):
        logger.addHandler(handler)
        _handlers.append(handler)
    logger.setLevel(logging.DEBUG)


def reset_logging() -> None:
    """Remove and close the handlers added by `setup_logging`."""
    for handler in _handlers:
        logger.removeHandler(handler)
        handler.close()
    _handlers.clear()
    logger.setLevel(logging.NOTSET)
