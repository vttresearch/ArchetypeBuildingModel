from setuptools import setup

from meta import __version__, description

setup(
    name="ArBuWe",
    python_requires=">=3.6",
    py_modules=["src/ArBuWe", "meta"],
    install_requires=[
        "geopandas",
        "atlite",
        "rioxarray",
        "numpy",
        "xarray",
        "matplotlib",
        "scipy",
        "pandas",
    ],
    version=__version__,
    description=description,
    author="Topi Rasku",
    author_email="topi.rasku@vtt.fi",
    license="",
)
