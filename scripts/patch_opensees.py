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
# 3b. sparsePython SOEs call find_package(Python3), which finds Termux's own
#     Python (3.14) and links/compiles against it. Use the injected target
#     Python headers instead and leave libpython to the host app.
for rel in ("SRC/system_of_eqn/eigenSOE/sparsePython/CMakeLists.txt",
            "SRC/system_of_eqn/linearSOE/sparsePython/CMakeLists.txt"):
    f = sys.argv[1] + "/" + rel
    t = open(f).read()
    t = t.replace("find_package(Python3 COMPONENTS Development REQUIRED)", "")
    t = t.replace("target_link_libraries(OPS_SysOfEqnPython PRIVATE Python3::Python)",
                  "target_include_directories(OPS_SysOfEqnPython PRIVATE ${Python_INCLUDE_DIRS})")
    open(f, "w").write(t); print("fixed", rel)
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
# OpenSeesReliabilityCommands.cpp: passes const char* to char* ctor params (clang error)
R = "SRC/interpreter/OpenSeesReliabilityCommands.cpp"
sub(R, "theGradientEvaluator, 0, filename);", "theGradientEvaluator, 0, (char*)filename);")
sub(R, "theProbabilityTransformation, filename,\n        relSensTag);", "theProbabilityTransformation, (char*)filename,\n        relSensTag);")
sub(R, "theFORMAnalysis, theFindCurvatures, filename);", "theFORMAnalysis, theFindCurvatures, (char*)filename);")
sub(R, "printFlag, filename, analysisTypeTag);", "printFlag, (char*)filename, analysisTypeTag);")
