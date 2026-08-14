"""Public Python API for XPID."""

from ._version import __version__
from .api import XPIDError, detect

__all__ = ["XPIDError", "detect", "__version__"]
