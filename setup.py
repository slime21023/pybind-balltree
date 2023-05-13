from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

__version__ = "0.0.1"

ext_modules = [
    Pybind11Extension(
        "pybind_balltree",
        sorted(glob("src/*.cpp")),  # Sort source files for reproducibility
    ),
]

setup(
    name="pybind_balltree",
    ext_modules=ext_modules,
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    cmdclass={"build_ext": build_ext},
)