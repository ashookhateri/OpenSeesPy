# Force a pure "py3-none-any" wheel: PyPI rejects linux_aarch64 tags, and
# Pydroid/Termux report a non-manylinux platform. The native modules are
# selected at import time by openseespyandroid/__init__.py.
from setuptools import setup
setup()
