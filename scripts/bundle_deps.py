"""Copy every non-Android-system shared library the modules need into libs/,
rewrite RUNPATH to $ORIGIN, and write a load order for the Python loader."""
import os, shutil, subprocess, sys
OUT = sys.argv[1]
PREFIX = os.environ.get("PREFIX", "/data/data/com.termux/files/usr")
LIBDIR = os.path.join(PREFIX, "lib")
SYSTEM = {"libc.so", "libm.so", "libdl.so", "liblog.so", "libz.so", "libandroid.so"}

def needed(path):
    out = subprocess.check_output(["readelf", "-d", path], text=True)
    return [l.split("[")[1].split("]")[0] for l in out.splitlines() if "(NEEDED)" in l]

seen, order = set(), []
def visit(lib):
    for n in needed(lib):
        if n in SYSTEM or n.startswith("libpython") or n in seen:
            continue
        src = os.path.join(LIBDIR, n)
        if not os.path.exists(src):
            print("WARNING: not found", n); continue
        seen.add(n)
        dst = os.path.join(OUT, "libs", n)
        shutil.copy2(os.path.realpath(src), dst)
        visit(dst)
        order.append(n)            # dependencies first

bins = [os.path.join(OUT, "bin", f) for f in os.listdir(os.path.join(OUT, "bin"))]
for b in bins:
    visit(b)
for f in [os.path.join(OUT, "libs", n) for n in order]:
    subprocess.check_call(["patchelf", "--set-rpath", "$ORIGIN", f])
for b in bins:
    subprocess.check_call(["patchelf", "--set-rpath", "$ORIGIN/../libs", b])
open(os.path.join(OUT, "libs", "load_order.txt"), "w").write("\n".join(order) + "\n")
print("bundled:", order)
