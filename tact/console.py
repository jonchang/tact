"""Logging and progress reporting shared by the command-line tools."""

import logging
import sys

import click

# Speed up logging for PyPy
logging._srcfile = None
logging.logThreads = False
logging.logProcesses = False
logging.logMultiprocessing = False

# Parent of the per-module loggers, so one set of handlers covers all of TACT
logger = logging.getLogger("tact")

_handlers: list[logging.Handler] = []
_console_level = logging.WARNING


def setup_logging(output: str, verbose: int) -> None:
    """Send TACT's log messages to `{output}.log.txt` and to standard error.

    The log file gets INFO messages, or DEBUG messages with `-vv`. Standard error gets warnings,
    INFO messages with `-v`, or DEBUG messages with `-vv`. Calling this again replaces the handlers
    from the previous call.

    Args:
        output: Base name of the output files.
        verbose: How many times `-v` was given.
    """
    global _console_level
    reset_logging()
    _console_level = logging.DEBUG if verbose >= 2 else logging.INFO if verbose == 1 else logging.WARNING
    file_handler = logging.FileHandler(output + ".log.txt")
    file_handler.setLevel(min(_console_level, logging.INFO))
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(_console_level)
    for handler in (file_handler, console_handler):
        logger.addHandler(handler)
        _handlers.append(handler)
    logger.setLevel(logging.DEBUG)


def reset_logging() -> None:
    """Remove and close the handlers added by `setup_logging`."""
    global _console_level
    for handler in _handlers:
        logger.removeHandler(handler)
        handler.close()
    _handlers.clear()
    logger.setLevel(logging.NOTSET)
    _console_level = logging.WARNING


def _isatty(stream) -> bool:
    try:
        return bool(stream.isatty())
    except (AttributeError, ValueError):
        return False


class Progress:
    """Progress through a fixed number of steps.

    Draws a progress bar on standard error when it is a terminal and log messages are not shown there.
    Otherwise it reports every tenth of the way, on standard error and in the log file, so that progress
    is still visible when the output is redirected or the bar would be broken up by log messages.

    Use `progressbar` to create one.
    """

    def __init__(self, iterable, label, length, item_show_func):
        """See `progressbar` for the arguments."""
        self.iterable = iterable
        self.label = label
        self.length = length
        self.pos = 0
        self._tenths = 0
        self._finished = False
        self._bar = None
        if _isatty(sys.stderr) and _console_level > logging.INFO:
            self._bar = click.progressbar(
                length=length,
                label=label,
                width=12,
                show_pos=True,
                item_show_func=item_show_func,
                file=sys.stderr,
            )
            self._bar.render_progress()

    def __enter__(self):
        """Return this object; the bar, if any, is already drawn."""
        return self

    def __exit__(self, *exc_info):
        """Finish the bar."""
        self.finish()

    def __iter__(self):
        """Yield each item of the iterable, advancing one step after each."""
        for item in self.iterable:
            yield item
            self.update(1, item)

    def update(self, n_steps: int, current_item=None) -> None:
        """Advance by `n_steps`, showing `current_item` next to the bar."""
        self.pos += n_steps
        if self._bar is not None:
            self._bar.update(n_steps, current_item)
        if self.length:
            tenths = min(self.pos * 10 // self.length, 10)
            if tenths > self._tenths:
                self._tenths = tenths
                self._report()

    def _report(self) -> None:
        message = f"{self.label}: {self.pos}/{self.length} ({self.pos * 100 // self.length}%)"
        logger.info(message)
        if self._bar is None and _console_level > logging.INFO:
            # the console handler does not show INFO messages, so print it there ourselves
            click.echo(message, err=True)

    def finish(self) -> None:
        """Finish drawing the bar, if there is one."""
        if self._bar is not None and not self._finished:
            self._bar.render_finish()
        self._finished = True


def progressbar(iterable=None, *, label: str, length: int, item_show_func=None) -> Progress:
    """Track progress over `iterable`, or over `length` steps advanced with `Progress.update`.

    Args:
        iterable: Items to iterate over, advancing one step per item.
        label: Name shown with the bar and at the start of each progress message.
        length: Total number of steps.
        item_show_func: Converts the current item to text shown next to the bar.

    Returns:
        A `Progress` that can be used as a context manager and iterated over.
    """
    return Progress(iterable, label, length, item_show_func)
