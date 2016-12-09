# adopted from http://peterdowns.com/posts/first-time-with-pypi.html
#from distutils.core import setup, find_packages
from setuptools import setup, find_packages

setup(
    name = 'cmap2d',
    packages = find_packages(exclude=['__pycache__']),  # this must be the same as the name above
    version = '0.0.1.dev0',
    description = 'A small library for generating and sampling 2D RGB colormaps.',
    author = 'Owen Derby',
    author_email = 'oderby@users.noreply.github.com',
    license = 'MIT',
    url = 'https://github.com/oderby/cmap2d',
    download_url = 'https://github.com/oderby/cmap2d/tarball/0.0.1.dev0',
    keywords = ['2D Colormaps', 'Colormaps', 'RGB Colormaps'],
    classifiers = [],
    install_requires = ['numpy', 'scipy', 'matplotlib'],
)
