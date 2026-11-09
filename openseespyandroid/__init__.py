"""openseespyandroid - OpenSeesPy for Android arm64 (Pydroid 3 / Termux).

Usage:
    import openseespyandroid.opensees as ops
    # or
    from openseespyandroid import opensees as ops
"""
import ctypes
import importlib.util
import os
import platform
import sys

__all__ = ["opensees"]

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIBS = os.path.join(_HERE, "libs")
_BIN = os.path.join(_HERE, "bin")


def _is_android():
    return (hasattr(sys, "getandroidapilevel")
            or "ANDROID_ROOT" in os.environ
            or os.path.exists("/system/build.prop"))


def _check_platform():
    if not _is_android():
        raise ImportError("openseespyandroid only runs on Android. "
                          "On desktop use: pip install openseespy")
    mach = platform.machine().lower()
    if mach not in ("aarch64", "arm64"):
        raise ImportError("openseespyandroid needs an arm64 (aarch64) device, "
                          "this one is %r" % mach)


def _preload_deps():
    # Load bundled shared libraries with RTLD_GLOBAL in dependency order so the
    # extension module finds them even where $ORIGIN RUNPATH is not honoured.
    order_file = os.path.join(_LIBS, "load_order.txt")
    if not os.path.exists(order_file):
        return
    with open(order_file) as f:
        names = [l.strip() for l in f if l.strip()]
    for name in names:
        path = os.path.join(_LIBS, name)
        if os.path.exists(path):
            ctypes.CDLL(path, mode=ctypes.RTLD_GLOBAL)


def _load():
    _check_platform()
    tag = "cp%d%d" % sys.version_info[:2]
    so = os.path.join(_BIN, "opensees_%s.so" % tag)
    if not os.path.exists(so):
        have = sorted(f[9:-3] for f in os.listdir(_BIN) if f.startswith("opensees_"))
        raise ImportError("No OpenSees build for Python %d.%d. Available: %s"
                          % (sys.version_info[0], sys.version_info[1], ", ".join(have)))
    _preload_deps()
    # The C module's init function is PyInit_opensees, so load it as "opensees".
    spec = importlib.util.spec_from_file_location("opensees", so)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    sys.modules["opensees"] = mod
    sys.modules[__name__ + ".opensees"] = mod
    return mod


opensees = _load()
