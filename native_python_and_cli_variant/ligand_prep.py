import os
import subprocess
import multiprocessing
import math
import shutil
import sys
import time
from rdkit import Chem
from rdkit.Chem import Descriptors

class LigandPrepper:
    def __init__(self, output_dir, cpu_count=1):
        self.output_dir = output_dir
        self.cpu_count = cpu_count
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if not shutil.which("obabel"):
            print("[CRITICAL ERROR] 'obabel' executable not found in PATH.")
            print("                 Please ensure OpenBabel is installed in the container.")
            sys.exit(1)

    def _split_structure_file(self, input_file, chunks_dir, num_chunks):
        """
        Splits SDF or MOL2 files into chunks.
        """
        if not os.path.exists(chunks_dir): os.makedirs(chunks_dir)
        
        filename = os.path.basename(input_file)
        ext = os.path.splitext(filename)[1].lower()
        
        if 'mol2' in ext:
            delimiter = "@<TRIPOS>MOLECULE"
            is_sdf = False
        else:
            delimiter = "$$$$"
            is_sdf = True

        molecules = []
        current_mol = []
        
        with open(input_file, 'r', errors='ignore') as f:
            for line in f:
                if delimiter in line:
                    if is_sdf:
                        current_mol.append(line)
                        molecules.append("".join(current_mol))
                        current_mol = []
                    else:
                        if current_mol: 
                            molecules.append("".join(current_mol))
                        current_mol = [line]
                else:
                    current_mol.append(line)
            if current_mol:
                molecules.append("".join(current_mol))

        total = len(molecules)
        if total == 0: return []
        if num_chunks > total: num_chunks = total
            
        chunk_size = math.ceil(total / num_chunks)
        chunk_files = []
        
        for i in range(num_chunks):
            start = i * chunk_size
            end = start + chunk_size
            subset = molecules[start:end]
            if not subset: continue
                
            fname = os.path.join(chunks_dir, f"chunk_{i}{ext}")
            with open(fname, 'w') as f:
                for m in subset: f.write(m)
                # Ensure newline to prevent hang
                f.write("\n")
            chunk_files.append(fname)
            
        return chunk_files

    def _run_obabel_chunk(self, args):
        """
        Worker that processes a chunk LIGAND-BY-LIGAND.
        This prevents one bad ligand from hanging the whole chunk for 30s.
        """
        input_chunk, output_chunk = args
        
        # 1. Read the chunk
        try:
            with open(input_chunk, 'r') as f:
                content = f.read()
        except:
            return False

        # 2. Parse into individual ligands
        is_mol2 = input_chunk.endswith('.mol2')
        delim = "@<TRIPOS>MOLECULE" if is_mol2 else "$$$$"
        
        if is_mol2:
            # Mol2 split: Delimiter is at START
            raw_mols = [delim + m for m in content.split(delim) if m.strip()]
        else:
            # SDF split: Delimiter is at END
            raw_mols = []
            for m in content.split("$$$$"):
                if m.strip():
                    raw_mols.append(m.strip() + "\n$$$$\n")

        valid_count = 0
        
        # 3. Process individually
        with open(output_chunk, 'w') as out_f:
            for i, mol_data in enumerate(raw_mols):
                # Write single temporary ligand
                temp_in = f"{input_chunk}_temp_{i}.{'mol2' if is_mol2 else 'sdf'}"
                temp_out = f"{input_chunk}_temp_out_{i}.sdf" # Always output SDF
                
                with open(temp_in, 'w') as t:
                    t.write(mol_data)
        
                # -isdf/-imol2 : Input format
                # -osdf : Output format (Forces SDF)
                # -O : Output filename
                input_flag = "-imol2" if is_mol2 else "-isdf"
                
                cmd = [
                    "obabel", input_flag, temp_in, 
                    "-osdf", "-O", temp_out,
                    "-p", "7.4", "--partialcharge", "gasteiger", "--gen3d"
                ]

                try:
                    subprocess.run(
                        cmd, 
                        check=True, 
                        stdout=subprocess.DEVNULL, 
                        stderr=subprocess.DEVNULL, 
                        timeout=3
                    )
                    
                    if os.path.exists(temp_out):
                        with open(temp_out, 'r') as res:
                            out_f.write(res.read())
                            valid_count += 1
                
                except subprocess.TimeoutExpired:
                    pass
                except Exception:
                    pass
                finally:
                    # Cleanup temp files immediately
                    if os.path.exists(temp_in): os.remove(temp_in)
                    if os.path.exists(temp_out): os.remove(temp_out)

        return valid_count > 0

    def _sanitize_sdf(self, sdf_path):
        """Ensures the final merged SDF is clean for Smina."""
        if not os.path.exists(sdf_path): return

        print(f"      > Sanitizing output for Smina compatibility...")
        try:
            with open(sdf_path, 'r') as f:
                content = f.read()
            
            raw_blocks = content.split("$$$$")
            clean_blocks = []
            
            for block in raw_blocks:
                if block.strip(): 
                    clean_mol = block.lstrip() + "\n$$$$\n"
                    clean_blocks.append(clean_mol)
            
            if clean_blocks:
                with open(sdf_path, 'w') as f:
                    f.write("".join(clean_blocks))
            else:
                 print("[WARNING] Sanitization resulted in empty file.")
                
        except Exception as e:
            print(f"[WARNING] Sanitization failed: {e}")

    def convert_and_clean(self, input_file):
        """Ligand Preparation Logic."""
        filename = os.path.basename(input_file)
        name, ext = os.path.splitext(filename)
        final_output = os.path.join(self.output_dir, f"{name}_prepared.sdf")
        
        print(f"[PREP] Processing {filename}...")

        # Parallelize if CPUs > 1
        if self.cpu_count > 1:
            chunks_dir = os.path.join(self.output_dir, "prep_chunks")
            chunk_files = self._split_structure_file(input_file, chunks_dir, self.cpu_count)
            
            if not chunk_files:
                print("[ERROR] Input file appears empty.")
                return None
                
            print(f"      > Parallelizing across {len(chunk_files)} workers...")
            
            tasks = []
            out_chunks = []
            for i, chunk in enumerate(chunk_files):
                out_name = os.path.join(chunks_dir, f"out_{i}.sdf")
                out_chunks.append(out_name)
                tasks.append((chunk, out_name))
            
            ctx = multiprocessing.get_context('spawn')
            
            with ctx.Pool(self.cpu_count) as pool:
                results = pool.map(self._run_obabel_chunk, tasks)
            
            print(" Done.")

            print("      > Waiting 5s for disk sync...")
            time.sleep(5)
                
            print(f"      > Merging results...")
            valid_mols = 0
            with open(final_output, 'w') as outfile:
                for fname in out_chunks:
                    if os.path.exists(fname):
                        with open(fname, 'r') as infile:
                            data = infile.read()
                            outfile.write(data)
                            if len(data) > 10: valid_mols += 1
            
            shutil.rmtree(chunks_dir, ignore_errors=True)
            
            if valid_mols == 0:
                print("[ERROR] No valid ligands were produced. Check OpenBabel installation.")
                return None

        else:
            # Serial Mode
            print(f"      > Running in serial mode (Single CPU)...")
            success = self._run_obabel_chunk((input_file, final_output))
            if not success:
                print("[ERROR] Ligand preparation failed (Serial Mode).")
                return None

        # Sanitize output
        self._sanitize_sdf(final_output)

        return final_output

    def filter_lipinski(self, input_sdf):
        output_filtered = input_sdf.replace(".sdf", "_lipinski.sdf")
        print(f"[PREP] Applying Lipinski's Rule of 5 Filter...")
        
        suppl = Chem.SDMolSupplier(input_sdf, sanitize=True)
        writer = Chem.SDWriter(output_filtered)
        
        passed = 0
        total = 0
        
        for mol in suppl:
            if mol is None: continue
            total += 1
            try:
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                hbd = Descriptors.NumHDonors(mol)
                hba = Descriptors.NumHAcceptors(mol)
                if (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10):
                    writer.write(mol)
                    passed += 1
            except:
                pass

        writer.close()
        print(f"[PREP] Lipinski Filter: {passed}/{total} molecules passed.")
        return output_filtered

def run_ligand_prep(input_file, output_dir, apply_lipinski=False, cpu_count=0):
    if cpu_count > 0:
        use_cpus = cpu_count
    else:
        avail_cpus = multiprocessing.cpu_count()
        use_cpus = max(1, int(avail_cpus * 0.75))
    
    prepper = LigandPrepper(output_dir, cpu_count=use_cpus)
    
    prepared_sdf = prepper.convert_and_clean(input_file)
    if not prepared_sdf: return None
    
    if apply_lipinski:
        final_sdf = prepper.filter_lipinski(prepared_sdf)
        return final_sdf
    
    return prepared_sdf
