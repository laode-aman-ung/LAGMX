import os
import sys
import subprocess
import glob
import re
import shutil
import tkinter as tk
import parmed as pmd

# Opsi tambahan untuk gmx mdrun (mis. "-ntmpi 1 -ntomp 16 -nb gpu"), diisi
# dari gmx_config.txt agar pemilihan sumber daya tidak di-hardcode.
MDRUN_OPTS = ""

def run_command(command):
    try:
        print(f"Running command: {command}")
        result = subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during execution: {e}")
        raise

def resolve_gaff_settings(atom_type):
    """Samakan versi GAFF yang dipakai tleap dan parmchk2 dengan -at pada antechamber"""
    at = (atom_type or "").strip().lower()
    if at == "gaff2":
        return "leaprc.gaff2", "gaff2"
    if at == "gaff":
        return "leaprc.gaff", "gaff"
    print(f"Warning: atom_type '{atom_type}' has no matching leaprc, falling back to gaff")
    return "leaprc.gaff", "gaff"

def detect_net_charge(mol_file):
    """Tentukan muatan bersih formal sebuah ligand dari file inputnya.

    Untuk PDB dibaca kolom muatan formal (kolom 79-80, mis. "N1+", "O1-");
    untuk MOL2 dijumlahkan muatan parsialnya lalu dibulatkan. Nilai ini harus
    diteruskan ke antechamber lewat -nc, karena tanpa itu antechamber
    mengasumsikan molekul netral dan menghasilkan muatan yang salah untuk
    ligand bermuatan tanpa pesan kesalahan apa pun.
    """
    ext = os.path.splitext(mol_file)[1].lower()
    total = 0.0
    if ext == '.mol2':
        in_atoms = False
        with open(mol_file) as f:
            for line in f:
                if line.startswith('@<TRIPOS>ATOM'):
                    in_atoms = True; continue
                if line.startswith('@<TRIPOS>'):
                    in_atoms = False
                if in_atoms:
                    parts = line.split()
                    if len(parts) >= 9:
                        try:
                            total += float(parts[8])
                        except ValueError:
                            pass
        return int(round(total))

    with open(mol_file) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            field = line[78:80].strip()
            if not field:
                continue
            sign = -1 if '-' in field else (1 if '+' in field else 0)
            digits = ''.join(c for c in field if c.isdigit())
            total += sign * (int(digits) if digits else 1)
    return int(round(total))

def build_ligand_ids(mol2_files):
    """Petakan setiap file ligand ke satu nama moleculetype yang unik.

    Nama residu pada file input tidak dijamin unik (mis. tiga salinan ligan
    yang sama semuanya bernama 'TZE'). Nama yang bertabrakan akan diberi
    akhiran urut sehingga jumlah ligan berapa pun tetap menghasilkan
    moleculetype, file posre, dan entri [ molecules ] yang berbeda-beda.
    """
    ids = {}
    used = set()
    for mol2_file in mol2_files:
        base = extract_ligand_id_from_mol2(mol2_file)
        if not base:
            base = os.path.splitext(os.path.basename(mol2_file))[0]
        base = base.strip()
        lig_id = base
        n = 1
        while lig_id in used:
            n += 1
            lig_id = f"{base}{n}"
        used.add(lig_id)
        ids[mol2_file] = lig_id
        if lig_id != base:
            print(f"Warning: duplicate ligand name '{base}', {mol2_file} renamed to '{lig_id}'")
    return ids

def extract_ligand_id_from_mol2(mol2_file):
    lig_id = None
    with open(mol2_file, 'r') as file:
        for line in file:
            if line.startswith('HETATM'):                
                parts = line.split()
                if len(parts) >= 4:
                    lig_id = parts[3]
                break
    return lig_id

def clean_itp_file(itp_file):
    """Remove the [ defaults ], [ system ], and [ molecules ] sections from the .itp file"""
    with open(itp_file, 'r') as file:
        lines = file.readlines()
    
    with open(itp_file, 'w') as file:
        skip = False
        for line in lines:
            if line.strip() in ['[ defaults ]', '[ system ]', '[ molecules ]']:
                skip = True
            elif skip and line.strip() == '':
                skip = False
            elif not skip:
                file.write(line.replace(".top", ".itp"))
        
def create_complex_pdb(receptor_new_files, mol2_files):
    complex_data = ""
    print(receptor_new_files)    
    try:
        for receptor_new_file in receptor_new_files:
            with open(receptor_new_file, 'r') as rec_new:
                rec_data = rec_new.read()
                complex_data += rec_data

        for mol2_file in mol2_files:
            lig_name = os.path.splitext(os.path.basename(mol2_file))[0]
            # Pakai koordinat tulisan parmed agar urutan dan nama atom persis
            # sama dengan topologi ligand; file input asli bisa memakai nama
            # atom yang berbeda sehingga grompp menolak/memberi peringatan.
            ligand_pdb = f"{lig_name}/{lig_name}_GMX.pdb"
            with open(ligand_pdb, 'r') as ligand_file:
                ligand_data = ligand_file.read()
                complex_data += ligand_data

        drop = ("TITLE", "REMARK", "MODEL", "ENDMDL", "END", "CONECT", "MASTER", "CRYST1")
        lines = complex_data.split('\n')
        filtered_lines = [line for line in lines if not line.startswith(drop)]

        with open("complex.pdb", 'w') as output_file:
            output_file.write('\n'.join(filtered_lines))
       
    except Exception as e:
        print(f"Error occurred during create_complex_pdb execution: {e}")
        raise

def create_complex_itp(mol2_files, lig_ids=None):
    if lig_ids is None:
        lig_ids = build_ligand_ids(mol2_files)

    atomtype_lines = []
    seen_atomtypes = set()
    moleculetype_itp = ""

    for mol2_file in mol2_files:
        lig_id = lig_ids[mol2_file]
        lig_name = os.path.splitext(mol2_file)[0]
        ligand_itp = f"{lig_name}/{lig_name}_GMX.itp"
        with open(ligand_itp, 'r') as f:
            flag_atomtypes = False
            flag_moleculetype = False
            lig_moleculetype_itp = ""
            orig_id = None

            for line in f:
                if '[ atomtypes ]' in line:
                    flag_atomtypes = True
                if not line.strip():
                    flag_atomtypes = False
                if flag_atomtypes:
                    # Dedup berdasarkan nama tipe atom, bukan pencocokan
                    # substring, supaya aman untuk jumlah ligand berapa pun.
                    key = line.split()[0] if line.split() else None
                    if line.strip().startswith((';', '[')):
                        if line not in atomtype_lines:
                            atomtype_lines.append(line)
                    elif key and key not in seen_atomtypes:
                        seen_atomtypes.add(key)
                        atomtype_lines.append(line)

                if '[ moleculetype ]' in line:
                    flag_moleculetype = True
                if flag_moleculetype:
                    lig_moleculetype_itp += line
                    # Baris pertama non-komentar setelah [ moleculetype ]
                    # adalah nama molekul bawaan parmed.
                    if orig_id is None and '[ moleculetype ]' not in line:
                        stripped = line.strip()
                        if stripped and not stripped.startswith(';'):
                            orig_id = stripped.split()[0]

            # Beri nama moleculetype yang unik
            if orig_id and orig_id != lig_id:
                lig_moleculetype_itp = lig_moleculetype_itp.replace(orig_id, lig_id, 1)

            # Restraint ligand harus berada di dalam blok moleculetype ligand
            # yang bersangkutan, bukan di level atas topol.top. Jika ditaruh di
            # level atas, semua include hanya menempel pada moleculetype terakhir
            # sehingga restraint salah untuk sistem multi-ligand.
            if lig_moleculetype_itp and not lig_moleculetype_itp.endswith("\n"):
                lig_moleculetype_itp += "\n"
            lig_moleculetype_itp += "\n; Ligand position restraint\n#ifdef POSRES\n"
            lig_moleculetype_itp += f'#include "{lig_id}-posre.itp"\n'
            lig_moleculetype_itp += "#endif\n\n"

            moleculetype_itp += lig_moleculetype_itp

    complex_itp = "".join(atomtype_lines) + "\n" + moleculetype_itp

    with open('ligand.itp', 'w') as complex_itp_file:
        complex_itp_file.write(complex_itp)

def split_topol2itp(topol_files):
    for topol_next in topol_files:
        atommoleculetype_lines = []
 
        write_atommoleculetype = False

        with open(topol_next, 'r') as topol_next_file:
            topol_next_lines = topol_next_file.readlines()

            for line in topol_next_lines:
                if "moleculetype" in line:
                    write_atommoleculetype = True
 
                elif "; Include water topology" in line:
                     write_atommoleculetype = False

                if write_atommoleculetype:
                    atommoleculetype_lines.append(line)

        # Simpan atommoleculetype dalam file dengan nama yang sesuai
        atommoleculetype_file_name = topol_next.replace(".top", ".itp")
        with open(atommoleculetype_file_name, 'w') as atommoleculetype_file:
            atommoleculetype_file.writelines(atommoleculetype_lines)

        # Hapus bagian atommoleculetype dan posre dari file topol_next
        with open(topol_next, 'w') as topol_next_file:
            filtered_lines = []
            remove_next = False

            for line in topol_next_lines:
                if any(keyword in line for keyword in ("moleculetype", "; Include Position restraint file")):
                    remove_next = True
                if "; Include water topology" in line:
                    remove_next = False

                if not remove_next:
                    filtered_lines.append(line)

            # Tambahkan baris yang diminta setelah baris kedua dari ; Include forcefield parameters
            index_ff_params = filtered_lines.index("; Include forcefield parameters\n")
            filtered_lines.insert(index_ff_params + 2, f"\n; Include chain topologies\n#include \"{atommoleculetype_file_name}\"\n")
            topol_next_file.writelines(filtered_lines)        

def create_topol_top(mol2_files, topol_files, lig_ids=None):
    if lig_ids is None:
        lig_ids = build_ligand_ids(mol2_files)
    print(topol_files)
    if os.path.exists("topol.top"):
        os.remove("topol.top")
    
    if topol_files:
        top1 = topol_files[0]  
        with open(top1, 'r') as topol_file:
            topol_lines = topol_file.readlines()

        for topol_next in topol_files[1:]:
            with open(topol_next, 'r') as topol_next:
                topol_next_lines = topol_next.readlines()
                new_topol1 = []  
                new_topol2 = []  
                
                include_chain = False
                molecules_section = False
                started_molecule_section = False

                for line in topol_next_lines:
                    if include_chain:
                        if line.startswith("#include"):
                            new_topol1.append(line)
                        else:
                            include_chain = False
                    elif line.strip() == "; Include chain topologies":
                        include_chain = True
                        
                    elif molecules_section:
                        if not started_molecule_section:
                            started_molecule_section = True
                        else:
                            new_topol2.append(line)
                    elif line.strip() == "[ molecules ]":
                        molecules_section = True
                    
                if new_topol1:
                    index = topol_lines.index("; Include chain topologies\n")
                    topol_lines = topol_lines[:index + 1] + new_topol1 + topol_lines[index + 1:]

                if new_topol2:
                    topol_lines.extend(new_topol2)
                    
        with open("topol.top", 'w') as topol_file:
            topol_file.writelines(topol_lines)

    lig_lines1 = "; Include ligand parameters\n#include \"ligand.itp\"\n\n"
    lig_lines2= "" 

    for file in mol2_files:
        lig_id = lig_ids[file]
        lig_lines2 += f'{lig_id}\t\t\t1\n'

    with open("topol.top", 'r') as topol_file:
        topol_lines = topol_file.readlines()

    index1 = None

    for i, line in enumerate(topol_lines):
        if line.strip() == "; Include forcefield parameters":
            index1 = i
            break  # Keluar dari loop setelah menemukan baris yang sesuai

    if index1 is not None:
        # Include restraint ligand sudah berada di dalam ligand.itp
        topol_lines.insert(index1 + 3, lig_lines1)

    index_ff_params = topol_lines.index("; Include forcefield parameters\n")
    merged_lines = topol_lines[index_ff_params:]
    merged_lines.extend([lig_lines2])

    with open("topol.top", 'w') as topol_file:
        topol_file.writelines(merged_lines)

def verify_ligand_charge(itp_file, expected, lig_name):
    """Pastikan muatan total topologi ligand sama dengan muatan bersih yang diminta."""
    total = 0.0
    in_atoms = False
    with open(itp_file) as f:
        for line in f:
            st = line.strip()
            if st.startswith('['):
                in_atoms = st.startswith('[ atoms ]')
                continue
            if in_atoms and st and not st.startswith(';'):
                parts = st.split()
                if len(parts) >= 7:
                    try:
                        total += float(parts[6])
                    except ValueError:
                        pass
    if abs(total - expected) > 0.05:
        print(f"Warning: {lig_name} topology sums to {total:+.3f} e but net charge "
              f"{expected:+d} was requested; check the charge method")
    else:
        print(f"{lig_name} net charge verified: {total:+.3f} e (expected {expected:+d})")
    return total

def gentop_gmx(mol_file, charge, atom_type, complex_dir, net_charge=None):
        print(f"Generating topology for {mol_file}")        
        lig_name = os.path.splitext(os.path.basename(mol_file))[0]
        file_extension = os.path.splitext(mol_file)[1]
        output_dir = os.path.join(".", lig_name)
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)
        shutil.copy(mol_file, output_dir)
        
        os.chdir(output_dir)
        
        # Construct the antechamber command
        if net_charge is None:
            net_charge = detect_net_charge(f"{lig_name}{file_extension}")
        print(f"Net charge for {lig_name}: {net_charge:+d}")
        if file_extension == '.mol2':
            print("Running antechamber for mol2...")
            run_command(f"antechamber -i {lig_name}.mol2 -fi mol2 -o {lig_name}_NEW.mol2 -fo mol2 -c {charge} -nc {net_charge} -s 2 -at {atom_type}")
        elif file_extension == '.pdb':            
            print("Running antechamber for pdb...")
            run_command(f"antechamber -i {lig_name}.pdb -fi pdb -o {lig_name}_NEW.mol2 -fo mol2 -c {charge} -nc {net_charge} -s 2 -at {atom_type}")            
        else:
            print(f"Unsupported file format: {file_extension}. Skipping {mol_file}.")
            return None, None

        # Generate GMX topology using parmchk2 and tleap
        leaprc_gaff, parmchk_gaff = resolve_gaff_settings(atom_type)
        print("Running parmchk2...")        
        run_command(f"parmchk2 -i {lig_name}_NEW.mol2 -f mol2 -o {lig_name}.frcmod -s {parmchk_gaff}")
        tleap_input = f"""
source leaprc.protein.ff14SB
source {leaprc_gaff}
mol = loadmol2 {lig_name}_NEW.mol2
loadamberparams {lig_name}.frcmod
saveamberparm mol {lig_name}.prmtop {lig_name}.inpcrd 
quit
        """
        tleap_input_file = os.path.join("tleap.in")
        with open(tleap_input_file, 'w') as f:
            f.write(tleap_input)        
        print("Running tleap...")
        run_command(f"tleap -f {tleap_input_file}")

        # Load the generated topology and coordinate files using parmed
        amber = pmd.load_file(f'{lig_name}.prmtop', f'{lig_name}.inpcrd')

        # Save the files in different formats
        print("Saving molecule and topology...")
        amber.save(f'{lig_name}.top')
        amber.save(f'{lig_name}_NEW.pdb')
        amber.save(f'{lig_name}.gro')

        # Copy the .top file to .itp
        shutil.copy(f'{lig_name}.top', f'{lig_name}_GMX.itp')
        shutil.copy(f'{lig_name}_NEW.pdb', f'{lig_name}_GMX.pdb')
        shutil.copy(f'{lig_name}.gro', f'{lig_name}_GMX.gro')

        # Clean the .itp file
        clean_itp_file(f'{lig_name}_GMX.itp')
        verify_ligand_charge(f'{lig_name}_GMX.itp', net_charge, lig_name)
        print("DONE")
        os.chdir("..")

CURATION_MARK = "REMARK   pyAutoGMX curated with pdbfixer"

def chains_in_pdb(pdb_file):
    """Kumpulkan chain ID yang muncul pada record ATOM sebuah file PDB."""
    chains = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith('ATOM'):
                cid = line[21]
                if cid not in chains:
                    chains.append(cid)
    return chains

def sequence_identity(receptor_file, seqres_records):
    """Bandingkan residu teramati dengan urutan SEQRES rujukan.

    Penomoran residu PDB diasumsikan mulai dari 1, sehingga residu bernomor n
    sesuai dengan entri ke-n pada SEQRES. Nilai yang dikembalikan adalah
    fraksi residu teramati yang namanya cocok.
    """
    seq = {}
    for chain in {l[11] for l in seqres_records}:
        names = []
        for l in seqres_records:
            if l[11] == chain:
                names.extend(l[19:].split())
        seq[chain] = names

    observed = {}
    with open(receptor_file) as f:
        for line in f:
            if line.startswith('ATOM'):
                observed[(line[21], int(line[22:26]))] = line[17:20].strip()

    hit = total = 0
    for (chain, num), name in observed.items():
        names = seq.get(chain)
        if not names or not (1 <= num <= len(names)):
            continue
        total += 1
        if names[num - 1] == name:
            hit += 1
    return hit / total if total else None

def attach_seqres(receptor_file, seqres_reference):
    """Sisipkan record SEQRES dari struktur rujukan ke file reseptor.

    pdbfixer hanya dapat memodelkan residu yang hilang bila urutan lengkapnya
    diketahui, dan urutan itu dibaca dari record SEQRES. Banyak file hasil
    ekspor program pemodelan membuang SEQRES, sehingga celah pada rantai tidak
    terdeteksi dan pdb2gmx akan menyambung dua residu yang berjauhan dengan
    ikatan peptida palsu. Fungsi ini mengambil SEQRES untuk chain ID yang benar
    dari struktur rujukan dan menempelkannya kembali.
    """
    with open(receptor_file) as f:
        content = f.read()
    if "SEQRES" in content:
        return False
    if not seqres_reference or not os.path.exists(seqres_reference):
        print(f"Warning: no SEQRES in {receptor_file} and no seqres_reference "
              f"available; missing residues cannot be modelled")
        return False

    wanted = set(chains_in_pdb(receptor_file))
    with open(seqres_reference) as f:
        records = [l for l in f if l.startswith("SEQRES") and l[11] in wanted]
    if not records:
        print(f"Warning: {seqres_reference} has no SEQRES for chains "
              f"{sorted(wanted)}")
        return False

    # Menempelkan urutan dari protein yang salah akan membuat pdbfixer
    # membangun residu yang keliru, jadi rujukan diverifikasi dulu terhadap
    # residu yang benar-benar teramati.
    ident = sequence_identity(receptor_file, records)
    if ident is None:
        print(f"Warning: cannot compare {receptor_file} with {seqres_reference}")
        return False
    if ident < 0.90:
        print(f"Warning: {os.path.basename(seqres_reference)} matches only "
              f"{ident:.0%} of the observed residues in "
              f"{os.path.basename(receptor_file)}; SEQRES not attached")
        return False
    print(f"SEQRES reference matches {ident:.0%} of observed residues")

    body = [l for l in content.splitlines(keepends=True)
            if not l.startswith(("SEQRES", "END"))]
    with open(receptor_file, 'w') as f:
        f.writelines(records)
        f.writelines(body)
        f.write("END\n")
    print(f"Attached {len(records)} SEQRES records for chains "
          f"{sorted(wanted)} to {receptor_file}")
    return True

def fix_receptor_structure(receptor_file, fixer_python, fixer_script,
                           seqres_reference=None):
    """Lengkapi residu internal yang hilang pada reseptor memakai pdbfixer.

    Dijalankan lewat interpreter terpisah (fixer_python) supaya pdbfixer/OpenMM
    tidak perlu dipasang di environment yang sama dengan parmed dan AmberTools.
    File yang sudah dikurasi ditandai dengan REMARK sehingga aman dijalankan
    ulang. Residu yang hilang di ujung rantai sengaja tidak dimodelkan agar
    rantai tidak diperpanjang melampaui daerah yang terpecahkan secara
    eksperimen.
    """
    with open(receptor_file) as f:
        head = f.read(4096)
    if CURATION_MARK in head:
        print(f"{receptor_file} already curated, skipping")
        return

    attach_seqres(receptor_file, seqres_reference)

    tmp = receptor_file + ".curated"
    run_command(f'"{fixer_python}" "{fixer_script}" "{receptor_file}" "{tmp}"')
    with open(tmp) as f:
        body = f.read()
    with open(receptor_file, 'w') as f:
        f.write(CURATION_MARK + "\n" + body)
    os.remove(tmp)

def run_pdb2gmx(receptor_file, forcefield, water):
    rec_name = os.path.splitext(receptor_file)[0]
    # -ff/-water membuat pdb2gmx berjalan tanpa prompt; -i memberi nama file
    # restraint per reseptor supaya beberapa rantai tidak saling menimpa
    # posre.itp yang sama.
    cmd = (f'gmx pdb2gmx -f {receptor_file} -o NEW_{rec_name}.pdb '
           f'-p topol_{rec_name}.top -i posre_{rec_name}.itp '
           f'-ff {forcefield} -water {water} -ignh')
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during pdb2gmx execution: {e}")
        raise

def generate_box_pdb(box_pdb, box_type, distance):
    try:
        subprocess.run(f'gmx editconf -f complex.pdb -o {box_pdb} -bt {box_type} -d {distance} -c', shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during generate_box_pdb execution: {e}")
        raise

def solvate_system(box_pdb, spc216_gro, topol_top, solv_gro):
    try:
        subprocess.run(f'gmx solvate -cp {box_pdb} -cs {spc216_gro} -p {topol_top} -o {solv_gro}', shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during solvate_system execution: {e}")
        raise

def ionization(solv_gro, topol_top, ions_gro):
    try:
        subprocess.run(f'gmx grompp -f ions.mdp -c {solv_gro} -p {topol_top} -o ions.tpr -maxwarn 1', shell=True, check=True)    
        subprocess.run(f'printf "SOL\\n" | gmx genion -s ions.tpr -o {ions_gro} -p {topol_top} -pname NA -nname CL -neutral', shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during ionization execution: {e}")
        raise

def minimization(ions_gro, topol_top):
    try:
        subprocess.run(f"gmx grompp -f em.mdp -c {ions_gro} -p {topol_top} -o em.tpr -maxwarn 1", shell=True, check=True)
        subprocess.run(f"gmx mdrun -v -deffnm em {MDRUN_OPTS}", shell=True, check=True)
        subprocess.run(f'echo "10\n0\n" | gmx energy -f em.edr -o potential.xvg', shell=True, check=True)    
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during minimization execution: {e}")
        raise
 
def make_index_system(em_gro, index_ndx, merged_group="Protein_LIG"):
    """Buat index.ndx berisi grup gabungan Protein + seluruh ligand.

    Grup bawaan "Other" berisi semua residu non-protein, non-air, non-ion,
    yaitu tepat seluruh ligand. Dengan memakai "Other" alih-alih menyebut
    nama tiap ligand, satu perintah yang sama berlaku untuk jumlah ligand
    berapa pun tanpa perlu tahu nama atau nomor grupnya.
    """
    try:
        # Pass 1: baca daftar grup bawaan untuk mengetahui nomor grup baru
        probe = subprocess.run(f'printf "q\\n" | gmx make_ndx -f {em_gro} -o {index_ndx}',
                               shell=True, check=True, capture_output=True, text=True)
        listing = probe.stdout + probe.stderr
        groups = re.findall(r'^\s*(\d+)\s+(\S+)\s*:', listing, re.MULTILINE)
        if not groups:
            raise RuntimeError("make_ndx did not report any default groups")
        names = [g[1] for g in groups]
        new_index = max(int(g[0]) for g in groups) + 1

        if "Other" not in names:
            raise RuntimeError("group 'Other' not found; no ligand detected in the system")

        # Pass 2: gabungkan Protein dengan seluruh ligand lalu beri nama
        select = f'"Protein" | "Other"\\nname {new_index} {merged_group}\\nq\\n'
        subprocess.run(f'printf \'{select}\' | gmx make_ndx -f {em_gro} -o {index_ndx}',
                       shell=True, check=True)
        print(f"Index group {new_index} created as {merged_group}")
    except (subprocess.CalledProcessError, RuntimeError) as e:
        print(f"Error occurred during make_index_system execution: {e}")
        raise

def replace_text_in_files(mol2_files, merged_group="Protein_LIG"):
    # Satu nama grup tetap dipakai berapa pun jumlah ligand, sehingga tc-grps
    # tidak lagi bergantung pada urutan penemuan file.
    mdp_files = sorted(glob.glob(os.path.join(".", "*.mdp")))
    for mdp_file in mdp_files:
        with open(mdp_file, 'r') as file:
            file_data = file.read()            
            file_data = file_data.replace("Protein_UNK", merged_group)
        with open(mdp_file, 'w') as file:
            file.write(file_data)
    
def process_ligand_restraint(mol2_files, lig_ids=None):
    if lig_ids is None:
        lig_ids = build_ligand_ids(mol2_files)
    for mol2_file in mol2_files:
        lig_name = os.path.splitext(os.path.basename(mol2_file))[0]
        lig_id = lig_ids[mol2_file]
        ligand_new_pdbs = glob.glob(f"{lig_name}/{lig_name}_GMX.pdb")
        
        for ligand_new_pdb in ligand_new_pdbs:   
            try:
                subprocess.run(f'printf "2 & ! a H*\\nq\\n" | gmx make_ndx -f {ligand_new_pdb} -o {lig_id}-index.ndx', shell=True, check=True)
                subprocess.run(f'printf "3\\n" | gmx genrestr -f {ligand_new_pdb} -n {lig_id}-index.ndx -o {lig_id}-posre.itp -fc 1000 1000 1000', shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error occurred during process_ligand_restraint execution: {e}")
                raise
                
def run_nvt_simulation():
    try:
        subprocess.run(f"gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 5", shell=True, check=True)
        subprocess.run(f"gmx mdrun -v -s nvt.tpr -deffnm nvt {MDRUN_OPTS}", shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during run_nvt_simulation execution: {e}")
        raise
    
def run_npt_simulation():
    try:
        subprocess.run(f"gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 5", shell=True, check=True)
        subprocess.run(f"gmx mdrun -v -s npt.tpr -deffnm npt {MDRUN_OPTS}", shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during run_npt_simulation execution: {e}")
        raise

def detect_stage(complex_dir_path):
    """Tentukan tahap terjauh yang SUDAH selesai untuk satu folder kompleks,
    supaya rerun bisa melanjutkan dari situ, bukan cuma dari tahap produksi.

    Hanya menghitung tahap yang file OUTPUT AKHIRNYA sudah ada (benar-benar
    selesai). Kalau run terputus di tengah suatu tahap (mis. NVT belum
    sempat menulis nvt.gro), tahap itu sendiri diulang dari checkpoint tahap
    sebelumnya -- tapi topologi/EM yang jauh lebih mahal tidak ikut diulang.
    """
    def exists(*names):
        return all(os.path.exists(os.path.join(complex_dir_path, n)) for n in names)
    if exists('md.tpr', 'md.cpt', 'npt.gro'):
        return 'production'
    if exists('npt.gro', 'npt.cpt', 'topol.top', 'index.ndx'):
        return 'npt_done'
    if exists('nvt.gro', 'nvt.cpt', 'topol.top', 'index.ndx'):
        return 'nvt_done'
    if exists('em.gro', 'topol.top', 'index.ndx'):
        return 'em_done'
    return 'none'

def run_production_simulation(target_ns=None, resume=False):
    try:
        if resume and os.path.exists("md.tpr") and os.path.exists("md.cpt"):
            # Hasil produksi sebelumnya sudah ada: jangan grompp ulang (itu akan
            # membangun tpr baru dari npt.gro/npt.cpt seolah mulai dari 0 ns lagi).
            # Cukup lanjutkan tpr yang ada.
            print(f"Melanjutkan produksi dari checkpoint menuju {target_ns} ns...")
        else:
            subprocess.run(f"gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr -maxwarn 5", shell=True, check=True)
        if target_ns is not None:
            # PENTING: pakai `convert-tpr -until` (waktu akhir ABSOLUT, ps), bukan
            # `mdrun -nsteps`/`convert-tpr -nsteps` yang keduanya berarti "step
            # TAMBAHAN dari checkpoint saat ini" -- kalau dipakai untuk target
            # absolut, hasilnya jadi checkpoint_step + override (run kelebihan
            # panjang), bukan sampai target_ns yang diminta.
            target_ps = target_ns * 1000
            subprocess.run(f"gmx convert-tpr -s md.tpr -until {target_ps} -o md.tpr", shell=True, check=True)
        # -cpi melanjutkan dari checkpoint bila run panjang terputus atau sengaja
        # diperpanjang (nsteps tpr, lewat -until di atas, > step yang sudah
        # tercapai); bila md.cpt belum ada, GROMACS memulai run baru seperti biasa.
        subprocess.run(f"gmx mdrun -v -s md.tpr -deffnm md -cpi md.cpt {MDRUN_OPTS}", shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during run_production_simulation execution: {e}")
        raise
        
if __name__ == "__main__":     
    current_directory = os.getcwd()
    with open('gmx_config.txt', 'r') as config_file:
            config_lines = config_file.readlines()
            config_variables = {}
            for line in config_lines:                
                    if not line.strip().startswith('#'):

                        try:
                            name, value = line.strip().split(': ')
                            name = name.strip()
                            value = value.strip()
                            config_variables[name] = value
                        except ValueError:
                            pass

    box_type = config_variables['box_type']
    box_type = box_type.lower()
    distance = config_variables['distance']
    solvent = config_variables['solvent']
    charge = config_variables['charge']
    atom_type = config_variables['atom_type']
    forcefield = config_variables.get('forcefield', 'amber03')
    water = config_variables.get('water', 'tip3p')
    merged_group = config_variables.get('merged_group', 'Protein_LIG')
    MDRUN_OPTS = config_variables.get('mdrun_options', '')
    fix_structure = config_variables.get('fix_structure', 'no').strip().lower() in ('yes', 'true', '1')
    fixer_python = config_variables.get('fixer_python', sys.executable)
    # fix_structure.py hidup satu tempat, di sebelah pyAutoGMX.py sendiri --
    # bukan diduplikasi ke tiap direktori kerja (current_directory bisa jadi
    # run_matrix/ atau folder lain, bukan tempat script ini berada).
    fixer_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fix_structure.py')
    seqres_reference = config_variables.get('seqres_reference', '').strip()
    if seqres_reference and not os.path.isabs(seqres_reference):
        seqres_reference = os.path.join(current_directory, seqres_reference)
    _nc = config_variables.get('net_charge', 'auto').strip().lower()
    net_charge_cfg = None if _nc in ('auto', '') else int(_nc)
    print(f"mdrun options: '{MDRUN_OPTS}'")

    # production_ns: panjang MD produksi (ns) yang diinginkan, menggantikan
    # edit manual nsteps di md.mdp. Bila kompleks sudah punya hasil produksi
    # sebelumnya (md.cpt), run itu DILANJUTKAN sampai production_ns tercapai,
    # bukan diulang dari 0 ns. Kosongkan/hapus key ini untuk pakai nsteps
    # bawaan md.mdp seperti sebelumnya (kompatibel dengan config lama).
    _production_ns_raw = config_variables.get('production_ns', '').strip()
    production_ns = float(_production_ns_raw) if _production_ns_raw else None
    if production_ns is not None:
        print(f"production target: {production_ns} ns ({production_ns * 1000} ps, absolut lewat convert-tpr -until)")
    else:
        print("production target: nsteps bawaan md.mdp (production_ns tidak diset)")

    mdps_files = sorted(glob.glob(os.path.join(current_directory, "*.mdp")))
    complex_dirs = sorted(d for d in os.listdir(current_directory)
                          if os.path.isdir(os.path.join(current_directory, d)) and d.startswith('complex'))
    print(f"Complex directories to process: {complex_dirs}")

    stage_map = {}
    resume_dirs = set()
    for complex_dir in complex_dirs:
        complex_dir_path = os.path.join(current_directory, complex_dir)
        stage = detect_stage(complex_dir_path) if production_ns is not None else 'none'
        stage_map[complex_dir] = stage
        if stage == 'production':
            # Hasil produksi sebelumnya ada: JANGAN dihapus, JANGAN diulang
            # dari system-prep/equilibrasi. Cukup lanjutkan produksinya nanti.
            resume_dirs.add(complex_dir)
            print(f"[{complex_dir}] hasil produksi sebelumnya terdeteksi -> akan dilanjutkan ke {production_ns} ns")
            continue
        if stage != 'none':
            # Prep/EM/NVT/NPT sebagian sudah ada -> jangan dihapus, lanjut
            # dari tahap itu (lihat loop di bawah).
            print(f"[{complex_dir}] tahap '{stage}' sebelumnya terdeteksi -> lanjut dari situ (tidak diulang dari 0)")
            continue
        for item in os.listdir(complex_dir_path):
            item_path = os.path.join(complex_dir_path, item)
            if item.startswith(('rec_', 'lig_', 'ref_')):
                continue

            if os.path.isdir(item_path):
                shutil.rmtree(item_path)
            else:
                os.remove(item_path)
        for mdp_file in mdps_files:
            shutil.copy(mdp_file, complex_dir_path)

    prepared_dirs = [d for d in complex_dirs if stage_map.get(d, 'none') != 'none']
    print("\n=== STAGE: system preparation for all complexes ===")
    for complex_dir in complex_dirs:
        if stage_map.get(complex_dir, 'none') != 'none':
            print(f"[{complex_dir}] skip system preparation (tahap '{stage_map[complex_dir]}' sudah ada)")
            continue
        os.chdir(complex_dir)
        # Kegagalan satu sistem tidak boleh menghentikan sistem lain
        try:
            # sorted() menjaga urutan reseptor dan ligand tetap sama antara
            # complex.pdb dan topol.top, dan membuat hasil reproducible.
            receptor_files = sorted(glob.glob('rec*.pdb'))
            mol2_files = sorted(glob.glob('lig*.pdb') + glob.glob('lig*.mol2'))
            lig_ids = build_ligand_ids(mol2_files)
            print(f"[{complex_dir}] receptors={receptor_files} ligands={lig_ids}")
            replace_text_in_files(mol2_files, merged_group)

            if receptor_files and mol2_files:
                if fix_structure:
                    local_refs = sorted(glob.glob('ref*.pdb'))
                    ref_here = local_refs[0] if local_refs else seqres_reference
                    for receptor_file in receptor_files:
                        fix_receptor_structure(receptor_file, fixer_python,
                                               fixer_script, ref_here)
                for mol2_file in mol2_files: 
                    gentop_gmx(mol2_file, charge, atom_type, complex_dir, net_charge_cfg)            
                for receptor_file in receptor_files:                
                    run_pdb2gmx(receptor_file, forcefield, water)
                topol_files = sorted(glob.glob('topol_*.top'))
                receptor_new_files = sorted(glob.glob('NEW*.pdb'))
                box_pdb = "box.pdb"
                topol_top = "topol.top"
                solv_gro = "solv.gro"
                ions_gro = "ions.gro"
                em_gro = "em.gro"
                index_ndx = "index.ndx"            
                create_complex_pdb(receptor_new_files, mol2_files)
                create_complex_itp(mol2_files, lig_ids)
                for topol_next in topol_files:
                    with open(topol_next, 'r') as topol_next_file:
                        topol_next_lines = topol_next_file.readlines()
                        if any("moleculetype" in line for line in topol_next_lines):
                            split_topol2itp([topol_next])
                create_topol_top(mol2_files, topol_files, lig_ids)
                generate_box_pdb(box_pdb, box_type, distance)
                solvate_system(box_pdb, solvent, topol_top, solv_gro)
                ionization(solv_gro, topol_top, ions_gro)
                minimization(ions_gro, topol_top)
                process_ligand_restraint(mol2_files, lig_ids)
                make_index_system(em_gro, index_ndx, merged_group)
                prepared_dirs.append(complex_dir)
        except Exception as e:
            print(f"[{complex_dir}] preparation failed: {e}")
        finally:
            os.chdir(current_directory)
       
    # Tahap ekuilibrasi dijalankan lebih dulu untuk SEMUA kompleks, baru
    # kemudian antrian produksi. Dengan begitu seluruh sistem sudah terekuilibrasi
    # sebelum satu pun run produksi yang panjang dimulai.
    equilibrated_dirs = []
    print("\n=== STAGE: NVT + NPT equilibration for all complexes ===")
    for complex_dir in prepared_dirs:
        stage = stage_map.get(complex_dir, 'none')
        if stage in ('npt_done', 'production'):
            print(f"[{complex_dir}] skip NVT+NPT (tahap '{stage}' sudah ada)")
            equilibrated_dirs.append(complex_dir)
            continue
        os.chdir(complex_dir)
        try:
            if stage == 'nvt_done':
                print(f"[{complex_dir}] skip NVT (sudah selesai sebelumnya)")
            else:
                print(f"\n--- [{complex_dir}] NVT ---")
                run_nvt_simulation()
            print(f"\n--- [{complex_dir}] NPT ---")
            run_npt_simulation()
            equilibrated_dirs.append(complex_dir)
        except Exception as e:
            print(f"[{complex_dir}] equilibration stage failed: {e}")
        finally:
            os.chdir(current_directory)

    print(f"\n=== STAGE: production MD queue {equilibrated_dirs} ===")
    for complex_dir in equilibrated_dirs:
        os.chdir(complex_dir)
        try:
            print(f"\n--- [{complex_dir}] production MD ---")
            run_production_simulation(target_ns=production_ns, resume=(complex_dir in resume_dirs))
        except Exception as e:
            print(f"[{complex_dir}] production stage failed: {e}")
        finally:
            os.chdir(current_directory)
