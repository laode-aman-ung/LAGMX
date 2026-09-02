# CLAUDE.md — LAGMX

Panduan kerja untuk sesi Claude Code di repositori ini.

---

## Ringkasan

LAGMX mengotomasi persiapan dan simulasi molecular dynamics GROMACS untuk
kompleks protein–ligan, tunggal maupun jamak. Menyiapkan satu sistem secara
manual butuh sembilan tahap, beberapa di antaranya interaktif; setiap ligan
tambahan melipatgandakan cara hal itu bisa salah — salinan ligan identik butuh
nama moleculetype berbeda atau GROMACS diam-diam merusak topologi, include
restraint harus berada di dalam blok moleculetype yang benar, tipe atom harus
dideduplikasi berdasarkan nama, dan muatan bersih tiap ligan harus dideteksi
lalu diteruskan dengan benar.

`LAGMX.py` menjalankan seluruhnya sebagai satu skrip dalam tiga tahap:
preparasi, ekuilibrasi (NVT lalu NPT), dan produksi.

**Sebelumnya bernama `pyAutoGMX`.** Diganti 2026-09-02 agar konsisten dengan
LADOCK dan LADEEP, sebelum pengajuan ulang ke JOSS. Tidak ada rujukan permanen
ke nama lama: belum pernah terbit di PyPI, belum ada DOI arsip. GitHub
mengalihkan URL lama, jadi tautan lama tetap hidup.

## Struktur

```
LAGMX.py           seluruh alat, satu skrip
fix_structure.py   pendamping: kurasi struktur reseptor
run_matrix/        satu-satunya contoh terpelihara dan bisa dijalankan:
                     gmx_config.txt + *.mdp + tujuh skenario complex_*/
                     run_lagmx.sh sebagai peluncur
tests/             pytest
paper.md/.bib      naskah JOSS
CITATION.cff       metadata sitasi
environment.yml    lingkungan conda bernama `lagmx`
```

## Hal terpenting sebelum mengubah apa pun

**Semua path dibaca relatif terhadap direktori tempat skrip dijalankan, bukan
tempat skripnya berada.** `gmx_config.txt`, berkas `.mdp`, dan direktori
`complex_*/` dicari di `os.getcwd()`. Menjalankan `python3 LAGMX.py` dari akar
repo akan gagal dengan `FileNotFoundError: gmx_config.txt` — itu perilaku
bawaan, bukan kerusakan.

Cara menjalankan yang benar:

```bash
cd run_matrix && ./run_lagmx.sh
```

Untuk proyek baru: salin `run_matrix/gmx_config.txt` dan `run_matrix/*.mdp` ke
direktori kerja Anda, taruh `complex_*/` di sebelahnya, lalu arahkan peluncur
ke `python3 ../LAGMX.py` dari dalam direktori itu.

## Cara menjalankan dan menguji

```bash
conda env create -f environment.yml && conda activate lagmx
```

```bash
python -m pytest tests/ -q
```

Dependensi non-Python ditarik conda: **GROMACS**, **AmberTools** (antechamber,
parmchk2), **ParmEd**. Semuanya harus ada di `PATH`; tidak ada yang dibundel.

## Konvensi yang teramati

- Satu skrip utama, bukan paket. Tidak ada `pyproject.toml`, tidak ada
  instalasi — dijalankan langsung dengan `python3`.
- `run_matrix/` adalah contoh **sekaligus** test matrix. Tujuh skenario
  menutupi ligan tunggal, tiga ligan, satu berkas gabungan, masukan mol2,
  ligan bermuatan, struktur rusak, dan multi-kompleks.
- Komentar kode dan pesan commit dalam bahasa Inggris.

## Yang perlu diperhatikan

1. **Naskah JOSS belum diajukan ulang.** Pengajuan sebagai `pyAutoGMX` gagal
   karena persoalan teknis di GitHub. Judul di `paper.md` dan `CITATION.cff`
   sudah memakai LAGMX.
2. **Tag lama v1.0.0–v1.0.3 dibuat sebelum penggantian nama** dan isinya masih
   menyebut pyAutoGMX. Tag tetap ada; rilis berikutnya sebaiknya versi baru
   yang menandai penggantian nama.
3. **Lisensi MIT** — berbeda dari LADOCK yang proprietary. Jangan menyamakan
   ketentuan kedua proyek.

## Data besar

Tidak ada. Repo ini 4,4 MB seluruhnya, dan `run_matrix/` hanya berisi kasus uji
kecil. Keluaran simulasi (`.xtc`, `.trr`, `.tpr`, `.edr`, `.cpt`) dihasilkan
saat dijalankan dan diabaikan `.gitignore` — tidak perlu Google Drive.

## Aturan sesi

### Awal sesi
1. Jalankan `git pull` sebelum melakukan apa pun.
2. Baca `STATE.md`.
3. Ringkas dalam 2–3 kalimat di mana pekerjaan terhenti.
4. Jangan mulai mengerjakan sebelum saya konfirmasi arahnya.

### Akhir sesi
Dipicu saat saya bilang "tutup sesi", "selesai", atau sejenisnya.
1. Perbarui `STATE.md`: perubahan, langkah berikutnya, yang macet. Ringkas,
   dengan tanggal dan nama mesin.
2. Perbarui `CLAUDE.md` bila ada keputusan atau konvensi baru.
3. Tampilkan daftar berkas yang akan di-commit dan tunggu persetujuan saya.
4. Commit dengan pesan deskriptif, lalu push.

### Sepanjang sesi
- Bahasa Indonesia untuk penjelasan, Inggris untuk komentar kode dan pesan
  commit.
- Jangan pernah menambahkan keluaran simulasi atau kredensial ke Git.
- Jangan membuat berkas baru bila mengedit yang ada sudah memadai.
