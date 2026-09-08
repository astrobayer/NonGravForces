#!/usr/bin/env bash
# runCChecks.sh -- everything about the C that can be checked without MATLAB.
#
# Run from the repository root:  bash mex/test/runCChecks.sh
#
# TWO CHECKS, AND WHY EACH EXISTS
#   1. The MEX gateway is compiled against mex/test/shim/mex.h. It does not link
#      and it is never run; the point is that a typo, a wrong argument count or
#      an undeclared identifier in ngfKernelGateway.c is found here rather than
#      on the owner's host, which is the only machine in this project that has a
#      MEX compiler.
#   2. The kernel is built and run against the frozen Python mirror fixtures in
#      tools/pyxcheck (WP-33), with the same bounds tests/tNgfMirror.m applies to
#      the MATLAB. That is the C's reference oracle, and it needs neither MATLAB
#      nor numpy, so the ubuntu audit job can gate it.
#
# Neither check replaces tests/tNgfMex.m, which measures the C against the MATLAB
# units directly and is the only place the 1e-14 parity claim can be made.
#
# Created 2026-09-08 by Claude Opus 5 (Anthropic).
set -euo pipefail

CC=${CC:-cc}
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

echo "runCChecks: compiler $($CC --version | head -1)"

echo "runCChecks: [1/2] compiling the MEX gateway against the shim (no link)"
$CC -std=c99 -Wall -Wextra -Werror \
    -I "$ROOT/mex/test/shim" -I "$ROOT/mex" \
    -c "$ROOT/mex/ngfKernelGateway.c" -o "$TMP/ngfKernelGateway.o"

echo "runCChecks: [2/2] building and running the parity harness"
$CC -std=c99 -O2 -Wall -Wextra -Werror \
    -I "$ROOT/mex" \
    -o "$TMP/ngfKernelParity" \
    "$ROOT/mex/test/ngfKernelParity.c" "$ROOT/mex/ngfKernel.c" -lm
"$TMP/ngfKernelParity" "$ROOT/tools/pyxcheck"

echo "runCChecks: clean"
