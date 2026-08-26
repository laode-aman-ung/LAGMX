"""Lengkapi residu/atom yang hilang pada file PDB reseptor dengan pdbfixer.

Hanya celah INTERNAL yang dimodelkan; residu yang hilang di ujung rantai
sengaja dibiarkan supaya rantai tidak diperpanjang melampaui daerah yang
benar-benar terpecahkan secara eksperimen. Hidrogen tidak ditambahkan karena
pdb2gmx yang akan menanganinya.
"""
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile

inp, out = sys.argv[1], sys.argv[2]
fixer = PDBFixer(filename=inp)

fixer.findMissingResidues()
chains = list(fixer.topology.chains())
internal = {}
for (chain_idx, res_idx), residues in fixer.missingResidues.items():
    chain = chains[chain_idx]
    n = len(list(chain.residues()))
    if res_idx == 0 or res_idx == n:
        print(f"  lewati {len(residues)} residu terminal pada rantai {chain.id}")
        continue
    internal[(chain_idx, res_idx)] = residues
    print(f"  modelkan {len(residues)} residu internal pada rantai {chain.id} "
          f"(setelah residu ke-{res_idx}): {' '.join(residues)}")
fixer.missingResidues = internal

fixer.findMissingAtoms()
n_at = sum(len(v) for v in fixer.missingAtoms.values())
n_te = sum(len(v) for v in fixer.missingTerminals.values())
if n_at or n_te:
    print(f"  lengkapi {n_at} atom berat hilang dan {n_te} atom terminal")
fixer.addMissingAtoms()

with open(out, 'w') as fh:
    PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
print(f"  -> {out}")
