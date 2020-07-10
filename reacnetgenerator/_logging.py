# cython: language_level=3
# cython: linetrace=True
"""Init logging."""
import logging
import coloredlogs
from ._version import __version__
from tqdm.autonotebook import tqdm

class TqdmLoggingHandler(logging.Handler):
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)

log = logging.getLogger(__name__)
log.addHandler(TqdmLoggingHandler())

coloredlogs.install(
    fmt=f'%(asctime)s - ReacNetGenerator {__version__} - %(levelname)s: %(message)s',
    loger=log,
    level=logging.INFO, milliseconds=True)

