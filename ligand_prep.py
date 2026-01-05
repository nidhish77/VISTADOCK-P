import os
import subprocess
import multiprocessing
import math
import shutil
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors

class LigandPrepper:
    def __init__(self, output_dir, cpu_count=1):
        self.output_dir = output_dir
        self.cpu_count = cpu_count
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def _count_mols(self, file_path):
        """Helper to count molecules in SDF or MOL2."""
        count = 0
        if not os.path.exists(file_path): return 0
        
        if file_path.endswith('.mol2'):
            delim = "@<TRIPOS>MOLECULE"
        else:
            delim = "$$$$" 
            
        try:
            with open(file_path, 'r', errors='ignore') as f:
                content = f.read()
                count = content.count(delim)
        except:
            return 0
        return count

    def _trim_to_last_delimiter(self, file_path):
        """
        Reads a file and truncates any incomplete data at the end.
        Ensures the file ends perfectly with '$$$$' or '@<TRIPOS>MOLECULE'.
        Prevents Smina 'tree.h(133)' crashes caused by killed processes.
        """
        if not os.path.exists(file_path): return
        
        is_mol2 = file_path.endswith('.mol2')
        delim = "@<TRIPOS>MOLECULE" if is_mol2 else "$$$$"
        
        with open(file_path, 'r', errors='ignore') as f:
            content = f.read()
            
        last_delim_pos = content.rfind(delim)
        
        if last_delim_pos == -1:
            with open(file_path, 'w') as f: f.write("") 
            return
        
        if not is_mol2:
            cutoff = last_delim_pos + 4
            clean_content = content[:cutoff] + "\n" 
        else:
            clean_content = content 

        with open(file_path, 'w') as f:
            f.write(clean_content)

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
                f.write("\n")
            chunk_files.append(fname)
            
        return chunk_files

    def _run_obabel_chunk(self, args):
        """
        Worker function with DYNAMIC TIMEOUT.
        """
        input_file, output_file = args
        
        input_count = self._count_mols(input_file)
        if input_count == 0: return False

        dynamic_timeout = (input_count * 30) + 60

        base_cmd = [
            "obabel", input_file, 
            "-osdf", "-O", output_file, 
            "-p", "7.4", "--partialcharge", "gasteiger"
        ]
        
        strategies = [
            ["--gen3d", "--minimize", "--steps", "500", "--ff", "MMFF94"],
            ["--gen3d", "--minimize", "--steps", "500", "--ff", "UFF"],
            [] 
        ]
        
        for flags in strategies:
            cmd = base_cmd + flags
            try:
                subprocess.run(
                    cmd, 
                    check=True, 
                    stdout=subprocess.DEVNULL, 
                    stderr=subprocess.DEVNULL, 
                    timeout=dynamic_timeout
                )
                
                # Success path
                if os.path.exists(output_file):
                    self._trim_to_last_delimiter(output_file) 
                    output_count = self._count_mols(output_file)
                    if output_count >= int(input_count * 0.9):
                        return True
                        
            except subprocess.TimeoutExpired:
                # Timeout path
                if os.path.exists(output_file):
                     self._trim_to_last_delimiter(output_file) 
                     output_count = self._count_mols(output_file)

                     if output_count >= int(input_count * 0.9):
                         return True
                continue
            except:
                continue
        return False

    def _sanitize_sdf(self, sdf_path):
        """
        Ensures the final merged SDF is clean for Smina.
        """
        if not os.path.exists(sdf_path): return

        print(f"      > Sanitizing output for Smina compatibility...")
        try:
            with open(sdf_path, 'r') as f:
                content = f.read()
            
            raw_blocks = content.split("$$$$")
            clean_blocks = []
            
            for block in raw_blocks:
                if block.strip(): 
                    # Strip leading newlines, add correct footer
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
        """
        Ligand Preparation Logic.
        """
        filename = os.path.basename(input_file)
        name, ext = os.path.splitext(filename)
        final_output = os.path.join(self.output_dir, f"{name}_prepared.sdf")
        
        print(f"[PREP] Processing {filename}...")

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
                for _ in pool.imap(self._run_obabel_chunk, tasks):
                    print(".", end="", flush=True)
            print(" Done.")
                
            print(f"      > Merging results...")
            with open(final_output, 'w') as outfile:
                for fname in out_chunks:
                    if os.path.exists(fname):
                        with open(fname, 'r') as infile:
                            outfile.write(infile.read())
            
            shutil.rmtree(chunks_dir, ignore_errors=True)

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