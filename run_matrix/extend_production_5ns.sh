#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Melanjutkan produksi MD dari 1 ns -> 5 ns untuk setiap complex_* di run_matrix.
#
# TIDAK memanggil pyAutoGMX.py (yang akan menghapus semua hasil 1 ns dan
# mengulang seluruh pipeline dari nol). Sebagai gantinya, script ini:
#   1. Memperpanjang md.tpr yang sudah ada lewat `gmx convert-tpr -extend`
#      (menambah 4000 ps -> total 5000 ps = 5 ns).
#   2. Melanjutkan mdrun dari checkpoint (md.cpt) dengan -cpi, sehingga
#      trajectory (md.xtc/md.edr/md.log) di-APPEND, bukan ditimpa.
#
# Jalankan HANYA setelah run_pyautogmx.sh selesai (semua complex_* punya
# md.gro + md.cpt hasil 1 ns).
#
# Pemakaian:
#   ./extend_production_5ns.sh              # semua complex_* yang sudah 1 ns
#   ./extend_production_5ns.sh complex_1single complex_2triple   # subset saja
# ---------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")"

# Uses whatever gmx is on PATH; activate the environment first, e.g.
#   conda activate pyautogmx
if ! command -v gmx >/dev/null 2>&1; then
    echo "GAGAL: 'gmx' tidak ada di PATH." >&2
    echo "Aktifkan environment-nya dulu: conda activate pyautogmx" >&2
    echo "(lihat ../environment.yml)" >&2
    exit 1
fi

export GMX_MAXBACKUP=-1
unset OMP_NUM_THREADS

TARGET_NS=5
EXTEND_PS=4000   # (5 - 1) ns * 1000 = 4000 ps ditambahkan ke tpr yang sudah 1 ns
MDRUN_OPTS="$(grep '^mdrun_options:' gmx_config.txt | cut -d: -f2- | sed 's/^ //')"

echo "=================== PREFLIGHT ==================="
echo "workdir      : $PWD"
echo "gmx          : $(command -v gmx)"
gmx --version 2>/dev/null | grep -E "GROMACS version|GPU support|SIMD instructions"
echo "mdrun opts   : $MDRUN_OPTS"
echo "target       : ${TARGET_NS} ns (extend +${EXTEND_PS} ps dari tpr yang ada)"
echo "-------------------------------------------------"
if ! gmx --version 2>/dev/null | grep -q "GPU support:.*CUDA"; then
    echo "GAGAL: gmx pada PATH tidak dibangun dengan dukungan CUDA." >&2
    exit 1
fi
if ! nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader; then
    echo "GAGAL: tidak ada GPU NVIDIA yang terdeteksi." >&2
    exit 1
fi
echo "================================================="
echo

targets=("$@")
if [ ${#targets[@]} -eq 0 ]; then
    mapfile -t targets < <(ls -d complex_*/ 2>/dev/null | sed 's#/##')
fi

overall_status=0
for d in "${targets[@]}"; do
    d=${d%/}
    echo "=== [$d] ==="
    if [ ! -f "$d/md.tpr" ] || [ ! -f "$d/md.cpt" ]; then
        echo "[$d] SKIP: md.tpr/md.cpt belum ada (produksi 1 ns belum selesai untuk kompleks ini)."
        overall_status=1
        continue
    fi

    ( cd "$d" && \
      cp -n md.tpr "md_1ns_backup.tpr" && \
      gmx convert-tpr -s md.tpr -extend "$EXTEND_PS" -o md.tpr && \
      gmx mdrun -v -s md.tpr -deffnm md -cpi md.cpt $MDRUN_OPTS
    ) 2>&1 | tee "$d/extend_5ns.log"
    status=${PIPESTATUS[0]}
    if [ "$status" -ne 0 ]; then
        echo "[$d] GAGAL (exit $status), lihat $d/extend_5ns.log"
        overall_status=1
    else
        echo "[$d] selesai diperpanjang ke ${TARGET_NS} ns."
    fi
    echo
done

echo "=================== RINGKASAN ==================="
printf '%-16s %-10s %s\n' SISTEM STATUS "durasi (ps) tercapai"
for d in "${targets[@]}"; do
    d=${d%/}
    ps_reached="-"
    if [ -f "$d/md.log" ]; then
        ps_reached=$(grep -oE '^ *Statistics over [0-9]+ steps.*' "$d/md.log" 2>/dev/null | tail -1)
    fi
    if [ -f "$d/md.gro" ] && [ -f "$d/md.tpr" ]; then
        printf '%-16s %-10s\n' "$d" "cek md.log"
    else
        printf '%-16s %-10s\n' "$d" "belum-selesai"
    fi
done
echo "==================================================="
exit $overall_status
