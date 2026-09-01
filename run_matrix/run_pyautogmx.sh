#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# pyAutoGMX launcher for the test matrix.
#
#   ./run_pyautogmx.sh                 # portable: use whatever gmx is on PATH
#   ./run_pyautogmx.sh --require-gpu   # refuse to run unless a CUDA GPU is used
#
# By default this runs anywhere: it uses the `gmx`, `antechamber` and `python3`
# already on your PATH (see ../environment.yml), and lets GROMACS pick the
# hardware it finds. Nothing about the simulation is set here -- every option
# is read from gmx_config.txt in this directory.
#
# --require-gpu reproduces the stricter behaviour used for the runs reported
# in the paper: the gmx binary must be built with CUDA and an NVIDIA GPU must
# be visible, so that mdrun fails loudly instead of silently falling back to
# the CPU. Combine it with an explicit `mdrun_options: ... -nb gpu` in
# gmx_config.txt to force GPU use for the non-bonded kernels.
# ---------------------------------------------------------------------------
set -uo pipefail
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$SELF")"

REQUIRE_GPU=0
for arg in "$@"; do
    case "$arg" in
        --require-gpu) REQUIRE_GPU=1 ;;
        -h|--help) sed -n '2,18p' "$SELF" | sed 's/^# \?//'; exit 0 ;;
        *) echo "Unknown option: $arg (try --help)" >&2; exit 2 ;;
    esac
done

export GMX_MAXBACKUP=-1
unset OMP_NUM_THREADS   # thread count comes from mdrun_options in gmx_config.txt

echo "=================== PREFLIGHT ==================="
echo "workdir     : $PWD"

missing=0
for tool in gmx antechamber tleap python3; do
    if command -v "$tool" >/dev/null 2>&1; then
        printf '%-12s: %s\n' "$tool" "$(command -v "$tool")"
    else
        printf '%-12s: MISSING\n' "$tool"; missing=1
    fi
done
if [ "$missing" -ne 0 ]; then
    cat >&2 <<'MSG'

Some requirements are missing. To install all of them:

    conda env create -f ../environment.yml
    conda activate pyautogmx

MSG
    exit 1
fi

gmx --version 2>/dev/null | grep -E "GROMACS version|GPU support|SIMD instructions"
echo "mdrun opts  : $(grep '^mdrun_options:' gmx_config.txt | cut -d: -f2-)"
echo "systems     : $(ls -d complex* 2>/dev/null | tr '\n' ' ')"

if [ "$REQUIRE_GPU" -eq 1 ]; then
    echo "-------------------------------------------------"
    if ! gmx --version 2>/dev/null | grep -q "GPU support:.*CUDA"; then
        echo "FAILED: the gmx on PATH is not built with CUDA support." >&2
        echo "Drop --require-gpu to run on the CPU instead." >&2
        exit 1
    fi
    if ! nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader; then
        echo "FAILED: no NVIDIA GPU detected." >&2
        echo "Drop --require-gpu to run on the CPU instead." >&2
        exit 1
    fi
    echo "GPU ready (--require-gpu)."
else
    if gmx --version 2>/dev/null | grep -q "GPU support:.*CUDA"; then
        echo "note        : CUDA build detected; GROMACS will use the GPU if it finds one."
    else
        echo "note        : CPU-only GROMACS build. This works, but production runs"
        echo "              will be slow. Pass --require-gpu to refuse to start"
        echo "              without a GPU."
    fi
fi
echo "================================================="
echo

python3 -u ../pyAutoGMX.py 2>&1 | tee run.log
status=${PIPESTATUS[0]}

echo
echo "=================== SUMMARY ====================="
echo "pyAutoGMX exit code : $status"
printf '%-16s %-9s %-9s %-9s %-9s %s\n' SYSTEM EM NVT NPT MD XTC
for d in complex*/; do
    printf '%-16s' "${d%/}"
    for f in em.gro nvt.gro npt.gro md.gro md.xtc; do
        if [ -f "$d$f" ]; then printf ' %-9s' "OK"; else printf ' %-9s' "-"; fi
    done
    echo
done
echo
if grep -qh -m1 -E "1 GPU selected|using .* GPU" complex*/md.log 2>/dev/null; then
    echo "--- GPU use per stage ---"
    grep -h -m1 -E "1 GPU selected|using .* GPU|PME tasks will do" complex*/md.log 2>/dev/null
fi
grep -h -E "^Performance:" complex*/md.log 2>/dev/null | sed 's/^/production -> /'
echo "================================================="
exit "$status"
