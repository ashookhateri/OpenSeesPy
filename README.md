# openseespyandroid

OpenSeesPy (OpenSees for Python) compiled for **Android arm64 (aarch64)**, for use in **Pydroid 3** and **Termux**.

## Install

```
pip install openseespyandroid
```

## Use

```python
import openseespyandroid.opensees as ops

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)
# ... same API as openseespy.opensees
```

Existing OpenSeesPy scripts: replace `import openseespy.opensees as ops` with
`import openseespyandroid.opensees as ops`.

## Supported

- Android 7.0+ (API 24), arm64-v8a devices
- Python 3.10, 3.11, 3.12, 3.13, 3.14 (the loader picks the matching build)
- Pydroid 3, Termux

## How it is built

`.github/workflows/build-release.yml` runs on a native GitHub arm64 runner inside
the official Termux Docker image (real Android bionic userland), compiles
OpenSees (C/C++ with clang, Fortran with LLVM flang), bundles the needed shared
libraries with `$ORIGIN` RUNPATH, and ships them in a `py3-none-any` wheel
(PyPI does not accept `linux_aarch64` wheels).

Pushing a tag `v*` builds, creates a GitHub Release and publishes to PyPI.

## License

OpenSees is (c) The Regents of the University of California. See `LICENSE`.
This package only redistributes compiled binaries of the upstream source.
