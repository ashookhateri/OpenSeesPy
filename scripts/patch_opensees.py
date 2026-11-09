"""Patch OpenSees CMakeLists.txt so it configures on Android/Termux (clang + flang)."""
import re, sys
p = sys.argv[1] + "/CMakeLists.txt"
s = open(p).read()
# 1. The "Clang" branch assumes macOS + gfortran static libs. Restrict it to APPLE.
s = s.replace('if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")',
              'if (APPLE AND "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")', 1)
# 2. No MPI / MKL on Android.
s = s.replace("find_package(MPI)", "# find_package(MPI) disabled for Android")
s = s.replace("find_package(MKL)", "# find_package(MKL) disabled for Android")
# 3. Let us inject the target Python (headers + stub libpython) from the command line.
open(p, "w").write(s)
print("patched", p)

# 4. Source fixes that LLVM flang rejects (gfortran only warns).
import os
root = sys.argv[1]
def sub(rel, old, new):
    f = os.path.join(root, rel)
    t = open(f).read()
    if old in t:
        open(f, "w").write(t.replace(old, new)); print("fixed", rel)
# SDM-UC.f: esspar declares oths(10) but reads oths(11); caller passes oths(12)
sub("SRC/material/nD/stressDensityModel/SDM-UC.f",
    "dimension props(nmats), oths(10)", "dimension props(nmats), oths(12)")
