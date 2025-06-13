############################################################
#
#   Init file for BSCAMPP, using the __init__.py from
#   SEPP as the original template. Current adaption comes
#   from https://github.com/c5shen/TIPP3.git
#
############################################################
from operator import itemgetter
import logging, os

# update system recursion limit to avoid issues
# not really needed for BSCAMPP but safe to update here
os.sys.setrecursionlimit(1000000)

__version__ = "1.0.8"
_INSTALL_PATH = __path__[0]

# global variables to store all loggers
__set_loggers = set()

# obtain the current logging level, default to INFO
def get_logging_level(logging_level='info'):
    logging_level_map = {
            'DEBUG': logging.DEBUG, 'INFO': logging.INFO,
            'WARNING': logging.WARNING, 'ERROR': logging.ERROR,
            'CRITICAL': logging.CRITICAL,
            }
    # obtain from environment variable to determine logging level, if
    # set by the user
    env_level = os.getenv('BSCAMPP_LOGGING_LEVEL')
    if env_level is not None:
        ll = env_level.upper()
    else:
        ll = logging_level.upper()
    # default to INFO if ll is not defined
    return logging_level_map.get(ll, logging.INFO)

# obtain a logger for a given file
def get_logger(name='bscampp', log_path=None, logging_level='info'):
    logger = logging.getLogger(name)
    if name not in __set_loggers:
        # set up a new logger for a name not in __set_loggers yet
        level = get_logging_level(logging_level)
        logging_formatter = logging.Formatter(
            ("[%(asctime)s] %(filename)s (line %(lineno)d):"
             " %(levelname) 8s: %(message)s"))
        logging_formatter.datefmt = "%H:%M:%S"
        logger.setLevel(level)

        # logging to stdout
        if log_path is None:
            ch = logging.StreamHandler()
        else:
            # use FileHandler for logging
            ch = logging.FileHandler(log_path, mode='a')
        ch.setLevel(level)
        ch.setFormatter(logging_formatter)
        logger.addHandler(ch)
        __set_loggers.add(name)
    return logger

# logging exception
def log_exception(logger):
    import traceback, io
    s = io.StringIO()
    traceback.print_exc(None, s)
    logger.error(s.getvalue())
    exit(1)
