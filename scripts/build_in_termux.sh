#!/data/data/com.termux/files/usr/bin/bash
# Builds opensees_cpXY.so for several Python versions inside a Termux aarch64
# environment (termux-docker on a GitHub arm64 runner, or a real phone).
set -euo pipefail
PY_VERSIONS="${PY_VERSIONS:-3.11.9 3.12.8 3.13.1}"
OPS_REF="${OPS_REF:-master}"
WORK="${WORK:-$HOME/aos}"
OUT="${OUT:-$WORK/out}"
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$WORK" "$OUT/bin" "$OUT/libs"

yes | pkg update
yes | pkg install clang flang cmake ninja make git python patchelf \
      libhdf5 zlib eigen tcl binutils curl

cd "$WORK"
[ -d OpenSees ] || git clone --depth 1 -b "$OPS_REF" https://github.com/OpenSees/OpenSees.git
python "$HERE/patch_opensees.py" "$WORK/OpenSees"

SYS_PYCONFIG="$(python -c 'import sysconfig;print(sysconfig.get_paths()["include"])')/pyconfig.h"

for V in $PY_VERSIONS; do
  MM="${V%.*}"; TAG="cp${MM/./}"
  echo "=== Python $MM ($TAG) ==="
  PYDIR="$WORK/py$MM"
  if [ ! -d "$PYDIR/Include" ]; then
    curl -L "https://www.python.org/ftp/python/$V/Python-$V.tgz" | tar xz -C "$WORK"
    mv "$WORK/Python-$V" "$PYDIR"
  fi
  # Headers: CPython Include/ + Termux's (Android aarch64) pyconfig.h
  cp "$SYS_PYCONFIG" "$PYDIR/Include/pyconfig.h"
  # Stub libpython with the soname Android Python builds use; real one is
  # provided by the host app (Pydroid/Termux) at load time.
  echo "" > "$PYDIR/stub.c"
  clang -shared -o "$PYDIR/libpython$MM.so" -Wl,-soname,"libpython$MM.so" "$PYDIR/stub.c"

  B="$WORK/build-$TAG"; rm -rf "$B"; mkdir -p "$B"; cd "$B"
  cmake "$WORK/OpenSees" -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_Fortran_COMPILER=flang \
    -DPython_FOUND=TRUE -DPython_INCLUDE_DIRS="$PYDIR/Include" \
    -DPython_LIBRARIES="$PYDIR/libpython$MM.so" \
    -DCMAKE_SHARED_LINKER_FLAGS="-Wl,--allow-shlib-undefined -Wl,-z,max-page-size=16384"
  ninja -j"$(nproc)" OpenSeesPy
  cp "$B/OpenSeesPy.so" "$OUT/bin/opensees_$TAG.so"
  cd "$WORK"
done

python "$HERE/bundle_deps.py" "$OUT"
echo "done -> $OUT"
