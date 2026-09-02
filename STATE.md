# Status

Diperbarui: 2026-09-02 — PC Ubuntu 24 (`arga`)

## Sedang dikerjakan

Penggantian nama `pyAutoGMX` → **LAGMX** selesai dan terdorong. Repo bersih.

## Yang berubah hari ini

- Nama diganti di seluruh proyek: 82 kemunculan pada 16 berkas, kedua ejaan
  (`pyAutoGMX` dan `pyautogmx`). Berkas ikut di-rename: `LAGMX.py`,
  `run_matrix/run_lagmx.sh`.
- Judul naskah JOSS dan `CITATION.cff` (title, repository-code, url) mengikuti.
- Repo GitHub di-rename ke `laode-aman-ung/LAGMX`; GitHub mengalihkan URL lama
  sehingga tautan yang sudah tersebar tetap hidup.
- Folder lokal pindah ke `~/riset/LAGMX`, mengikuti pola path identik lintas
  mesin.
- Ditambahkan berkas wajib sinkronisasi: `CLAUDE.md`, `STATE.md`,
  `.gitattributes`.

## Langkah berikutnya

1. **Ajukan ulang ke JOSS** dengan nama LAGMX. Pengajuan sebelumnya sebagai
   pyAutoGMX gagal karena persoalan teknis di GitHub.
2. Buat rilis baru yang menandai penggantian nama, lalu selaraskan `version`
   dan `date-released` di `CITATION.cff`. Saat ini masih `1.0.3` /
   `2026-09-01`, yaitu rilis terakhir dengan nama lama.
3. JOSS mensyaratkan arsip ber-DOI (Zenodo) saat penerimaan. Buat arsipnya
   **setelah** nama final, supaya DOI-nya mencatat LAGMX.

## Tertunda / macet

- Tag `v1.0.0`–`v1.0.3` dibuat sebelum penggantian nama dan isinya masih
  menyebut pyAutoGMX. Dibiarkan apa adanya sebagai catatan sejarah.
- Belum ada rilis dengan nama baru.

## Keputusan terakhir

- 2026-09-02 — Nama diganti **sebelum** pengajuan ulang JOSS, karena setelah
  makalah terbit nama akan terkunci di DOI dan di sitasi orang lain. Saat
  penggantian dilakukan, belum ada rujukan permanen: nihil di PyPI, nihil DOI
  arsip.
- 2026-09-02 — Tidak ada data besar dan tidak perlu Google Drive; keluaran
  simulasi dihasilkan saat dijalankan dan sudah diabaikan `.gitignore`.
- 2026-09-02 — Lisensi tetap **MIT**, berbeda dari LADOCK yang proprietary.
