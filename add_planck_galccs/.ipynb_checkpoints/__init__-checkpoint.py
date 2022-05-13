from .galcc_utils import *
from .generate_galcc_maps import *

# Enforce Python version check during package import.
import sys

__minimum_python_version__ = "3.6"

class UnsupportedPythonError(Exception):
    pass

if sys.version_info < tuple(int(val for val in __minimum_python_version__.split("."))):
    raise UnsupportedPythonError("Package does not support Python < {}".format(__minimum_python_version__))
    
try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    from importlib_metadata import version, PackageNotFoundError
    
try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "unknown"