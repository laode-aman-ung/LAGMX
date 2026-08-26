#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Launcher pyAutoGMX
#   GROMACS  : 2026.0 conda-forge (CUDA + AVX2_256) dari env "gromacs_env"
#   AmberTools / parmed : miniconda base
#   GPU      : dipaksa lewat "mdrun_options: ... -nb gpu" di gmx_config.txt,
#              sehingga mdrun GAGAL (bukan diam-diam pindah ke CPU) bila GPU
#              tidak terpakai.
# Semua opsi simulasi dibaca dari gmx_config.txt, tidak ada yang di-hardcode.
# ---------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")"

export PATH="$HOME/miniconda3/envs/gromacs_env/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/gromacs_env/lib:${LD_LIBRARY_PATH:-}"
export GMX_MAXBACKUP=-1
unset OMP_NUM_THREADS   # jumlah thread diatur lewat mdrun_options di gmx_config.txt

echo "=================== PREFLIGHT ==================="
echo "workdir     : $PWD"
echo "gmx         : $(command -v gmx)"
gmx --version 2>/dev/null | grep -E "GROMACS version|GPU support|SIMD instructions"
echo "antechamber : $(command -v antechamber)"
echo "tleap       : $(command -v tleap)"
echo "mdrun opts  : $(grep '^mdrun_options:' gmx_config.txt | cut -d: -f2-)"
echo "sistem      : $(ls -d complex* | tr '\n' ' ')"
echo "-------------------------------------------------"

# 1. binary gmx harus punya dukungan CUDA
if ! gmx --version 2>/dev/null | grep -q "GPU support:.*CUDA"; then
    echo "GAGAL: gmx pada PATH tidak dibangun dengan dukungan CUDA." >&2
    exit 1
fi
# 2. GPU harus benar-benar terlihat driver
if ! nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader; then
    echo "GAGAL: tidak ada GPU NVIDIA yang terdeteksi." >&2
    exit 1
fi
echo "GPU siap dipakai (mdrun dipaksa -nb gpu)."
echo "================================================="
echo

"$HOME/miniconda3/bin/python3" -u ../pyAutoGMX.py 2>&1 | tee run.log
status=${PIPESTATUS[0]}

echo
echo "=================== RINGKASAN ==================="
echo "exit code pyAutoGMX : $status"
printf '%-16s %-9s %-9s %-9s %-9s %s\n' SISTEM EM NVT NPT MD XTC
for d in complex*/; do
    printf '%-16s' "${d%/}"
    for f in em.gro nvt.gro npt.gro md.gro md.xtc; do
        if [ -f "$d$f" ]; then printf ' %-9s' "OK"; else printf ' %-9s' "-"; fi
    done
    echo
done
echo
echo "--- konfirmasi pemakaian GPU pada tiap tahap ---"
grep -h -m1 -E "1 GPU selected|using .* GPU|PME tasks will do" complex*/md.log 2>/dev/null
grep -h -E "^Performance:" complex*/md.log 2>/dev/null | sed 's/^/produksi 1 ns -> /'
echo "================================================="
