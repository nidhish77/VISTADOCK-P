import os
import json
import subprocess
import multiprocessing
import math
import shutil
from openff.toolkit.topology import Molecule
from mmgbsa_scorer import GBSAScorer
from plip_analyzer import PLIPAnalyzer
from openmm.app import Modeller, PDBFile 

# TERMINAL COLORS 
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def log(tag, message, color=Colors.ENDC):
    print(f"{color}[{tag.upper()}] {message}{Colors.ENDC}")

def ensure_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)

# --- HELPER: Count Molecules ---
def count_molecules(sdf_path):
    if not os.path.exists(sdf_path): return 0
    try:
        with open(sdf_path, 'r', errors='ignore') as f: 
            return f.read().count("$$$$")
    except:
        return 0

# --- HELPER: Split SDF (Robust Version) ---
def split_sdf_file(input_sdf, chunks_dir, num_chunks):
    """
    Splits a large SDF file into 'num_chunks' smaller files.
    [CRITICAL FIX]: Ensures no leading newlines in chunks, which cause Smina to crash.
    """
    if not os.path.exists(chunks_dir):
        os.makedirs(chunks_dir)
    
    # Read entire file content
    with open(input_sdf, 'r') as f:
        content = f.read()

    # 1. Split by the SDF delimiter '$$$$'
    raw_blocks = content.split("$$$$")
    
    clean_mols = []
    for block in raw_blocks:
        if block.strip(): 
            clean_mol = block.lstrip() + "\n$$$$\n"
            clean_mols.append(clean_mol)

    total_mols = len(clean_mols)
    if total_mols == 0: return []
    
    if num_chunks > total_mols:
        num_chunks = total_mols
        
    chunk_size = math.ceil(total_mols / num_chunks)
    chunk_files = []
    
    for i in range(num_chunks):
        start = i * chunk_size
        end = start + chunk_size
        subset = clean_mols[start:end]
        
        if not subset: continue
            
        chunk_name = os.path.join(chunks_dir, f"chunk_{i}.sdf")
        with open(chunk_name, 'w') as f:
            for mol_block in subset:
                f.write(mol_block)
        
        chunk_files.append(chunk_name)
        
    return chunk_files

def merge_sdfs(chunk_files, output_file):
    """Merges SDF chunks efficiently."""
    with open(output_file, 'w') as outfile:
        for fname in chunk_files:
            if os.path.exists(fname):
                with open(fname, 'r') as infile:
                    outfile.write(infile.read())

# --- RUNNER: Single Smina Command ---
def run_smina_cmd(receptor, ligand, output, exhaustiveness, scoring_func, cpu_count, center, size):
    cmd = [
        'smina', '-r', receptor, '-l', ligand, '-o', output,
        '--exhaustiveness', str(exhaustiveness), '--scoring', scoring_func,
        '--cpu', str(cpu_count),
        '--center_x', str(center[0]), '--center_y', str(center[1]), '--center_z', str(center[2]),
        '--size_x', str(size[0]), '--size_y', str(size[1]), '--size_z', str(size[2]),
        '--quiet' 
    ]
    try:
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError:
        return False

# --- RUNNER: Parallel Engine ---
def run_parallel_smina(receptor, input_sdf, output_sdf, exhaustiveness, scoring_func, total_cpus, center, size, work_dir):
    """
    Adaptive Execution Engine:
    - Small Batches (<50): Internal Threading (Low Overhead)
    - Large Batches (>50): Parallel Processes (High Throughput)
    """
    num_mols = count_molecules(input_sdf)
    
    if num_mols < 50 or total_cpus <= 1:
        log(scoring_func, f"   > Mode: INTERNAL THREADING ({num_mols} ligands, {total_cpus} CPUs)", Colors.BLUE)
        run_smina_cmd(receptor, input_sdf, output_sdf, exhaustiveness, scoring_func, total_cpus, center, size)
        return

    # --- Start Parallel Execution ---
    log(scoring_func, f"   > Mode: PARALLEL SPLIT ({num_mols} ligands over {total_cpus} workers)", Colors.BLUE)
    
    chunks_dir = os.path.join(work_dir, "chunks_in")
    out_dir = os.path.join(work_dir, "chunks_out")
    if not os.path.exists(out_dir): os.makedirs(out_dir)

    # 1. Split Inputs
    input_chunks = split_sdf_file(input_sdf, chunks_dir, total_cpus)
    
    if not input_chunks: return 

    processes = []
    output_chunks = []

    # 2. Launch Swarm
    for i, chunk_in in enumerate(input_chunks):
        chunk_out = os.path.join(out_dir, f"out_{i}.sdf")
        output_chunks.append(chunk_out)
        
        cmd = [
            'smina', '-r', receptor, '-l', chunk_in, '-o', chunk_out,
            '--exhaustiveness', str(exhaustiveness), '--scoring', scoring_func,
            '--cpu', '1',  # Force 1 CPU per worker
            '--center_x', str(center[0]), '--center_y', str(center[1]), '--center_z', str(center[2]),
            '--size_x', str(size[0]), '--size_y', str(size[1]), '--size_z', str(size[2]),
            '--quiet'
        ]
        p = subprocess.Popen(cmd)
        processes.append(p)

    # 3. Wait for completion
    for p in processes:
        p.wait()

    # 4. Merge
    merge_sdfs(output_chunks, output_sdf)
    
    # 5. Cleanup
    shutil.rmtree(chunks_dir, ignore_errors=True)
    shutil.rmtree(out_dir, ignore_errors=True)

# --- WORKER: Complex Generation (Single Item) ---
def _generate_complex_worker(data):
    """Helper for pool processing"""
    rec_path, mol, out_path = data
    try:
        receptor = PDBFile(rec_path)
        modeller = Modeller(receptor.topology, receptor.positions)
        lig_top = mol.to_topology().to_openmm()
        lig_pos = mol.conformers[0].to_openmm()
        modeller.add(lig_top, lig_pos)
        with open(out_path, 'w') as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f)
        return True
    except:
        return False

# --- MANAGER: Complex Generation (Parallelized) ---
def generate_complexes(rec_pdb_path, lig_sdf_path, output_dir, num_workers=1):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    try:
        # Load all molecules
        mols = Molecule.from_file(lig_sdf_path, allow_undefined_stereo=True)
        if isinstance(mols, Molecule): mols = [mols] 
        
        tasks = []
        for i, mol in enumerate(mols):
            raw_name = mol.name if mol.name else f"Ligand_{i+1}"
            safe_name = "".join([c for c in raw_name if c.isalnum() or c in ('_', '-')])
            complex_filename = os.path.join(output_dir, f"{safe_name}_complex.pdb")
            tasks.append((rec_pdb_path, mol, complex_filename))

        active_workers = min(num_workers, len(tasks))
        if active_workers < 1: active_workers = 1
        
        with multiprocessing.Pool(active_workers) as pool:
            results = pool.map(_generate_complex_worker, tasks)
            
        return sum(results)
    except Exception as e:
        print(f"[ERROR] Failed to generate complexes: {e}")
        return 0

def convert_pdb_to_pdbqt(pdb_path, output_pdbqt_path):
    cmd = ["obabel", "-ipdb", pdb_path, "-opdbqt", "-O", output_pdbqt_path, "-h", "-xr", "--partialcharge", "gasteiger"]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        return False
    except FileNotFoundError:
        print(f"{Colors.FAIL}[ERROR] 'obabel' not found.{Colors.ENDC}")
        return False

def filter_sdf(input_sdf, output_sdf, fraction):
    if not os.path.exists(input_sdf): return 0
    with open(input_sdf, 'r') as f: content = f.read()
    blocks = content.split("$$$$")
    best_unique_candidates = {}
    for block in blocks:
        if not block.strip(): continue
        lines = block.strip().split('\n')
        lig_id = lines[0].strip()
        score = 100.0
        for i, line in enumerate(lines):
            if "> <minimizedAffinity>" in line:
                try: score = float(lines[i+1].strip())
                except: pass
                break
        if lig_id not in best_unique_candidates or score < best_unique_candidates[lig_id]['score']:
            best_unique_candidates[lig_id] = {'score': score, 'content': block.strip() + '\n\n$$$$\n'}
    unique_data = sorted(list(best_unique_candidates.values()), key=lambda x: x['score'])
    keep_count = max(1, int(len(unique_data) * fraction))
    survivors = unique_data[:keep_count]
    with open(output_sdf, 'w') as f:
        for item in survivors: f.write(item['content'])
    return len(survivors)

def analyze_track_results(precise_survivors_path, rec_pdb_path, output_dir, tag_color, tag_name, run_gbsa=True, run_plip=True):
    results = {}
    gb_scorer = GBSAScorer(rec_pdb_path) if run_gbsa else None
    plip_tool = PLIPAnalyzer(rec_pdb_path) if run_plip else None
    plip_out_dir = os.path.join(output_dir, "plip_reports")
    
    try:
        mols = Molecule.from_file(precise_survivors_path, allow_undefined_stereo=True)
        if isinstance(mols, Molecule):
            mols = [mols]
    except:
        log(tag_name, "Error loading survivors for analysis.", Colors.FAIL)
        return {}

    msg = f"Analyzing {len(mols)} candidates..." if run_gbsa else f"Profiling interactions for {len(mols)} candidates (GBSA Skipped)..."
    log(tag_name, msg, tag_color)

    for i, mol in enumerate(mols):
        name = mol.name if mol.name else f"Ligand_{i+1}"
        smina_score = mol.properties.get('minimizedAffinity', "N/A")
        
        if run_gbsa:
            try: gbsa = gb_scorer.calculate_single_dG(mol)
            except Exception as e: 
                print(f"[{tag_name}] GBSA Error for {name}: {e}") 
                gbsa = -999.9 
        else:
            gbsa = "Skipped"
            
        if run_plip:
            try:
                interactions = plip_tool.analyze_interactions(mol, output_dir=plip_out_dir)
            except:
                interactions = "Error"
        else:
            interactions = "Skipped"
            
        results[name] = {'smina': float(smina_score) if smina_score != "N/A" else 0.0, 'gbsa': gbsa, 'plip': interactions}
    return results

def run_screening_track(scoring_func, base_output_dir, return_dict, config):
    if scoring_func == 'vina': color = Colors.GREEN
    else: color = Colors.CYAN
        
    assigned_cpus = config.get('cpu_count', 1)
    
    rec_pdb   = config['receptor_pdb']
    lig_sdf   = config['ligand']
    center = (config['center_x'], config['center_y'], config['center_z'])
    size   = (config['size_x'], config['size_y'], config['size_z'])

    track_dir = os.path.join(base_output_dir, f"{scoring_func}_track")
    ensure_directory(track_dir)
    log(scoring_func, f"Track Initiated using {assigned_cpus} CPUs in Data-Parallel Mode.", color)

    # 0. Prep Receptor
    rec_pdbqt = os.path.join(track_dir, "receptor_prepared.pdbqt")
    log(scoring_func, "Converting Receptor PDB to PDBQT...", color)
    if not convert_pdb_to_pdbqt(rec_pdb, rec_pdbqt):
        log(scoring_func, "Failed to convert Receptor. Aborting.", Colors.FAIL); return

    # 1. RAPID
    rapid_out = os.path.join(track_dir, "results_RAPID.sdf")
    rapid_surv = os.path.join(track_dir, "survivors_RAPID.sdf")
    exh_rapid = config.get('exh_rapid', 2); frac_rapid = config.get('frac_rapid', 0.5)
    
    log(scoring_func, f"Step 1: RAPID (Exh={exh_rapid}, Keep={frac_rapid*100}%)", color)
    run_parallel_smina(rec_pdbqt, lig_sdf, rapid_out, exh_rapid, scoring_func, assigned_cpus, center, size, track_dir)
    count = filter_sdf(rapid_out, rapid_surv, frac_rapid)
    if count == 0: return

    # 2. STANDARD
    standard_out = os.path.join(track_dir, "results_STANDARD.sdf")
    standard_surv = os.path.join(track_dir, "survivors_STANDARD.sdf")
    exh_standard = config.get('exh_standard', 8); frac_standard = config.get('frac_standard', 0.3)
    
    log(scoring_func, f"Step 2: STANDARD (Exh={exh_standard}, Keep={frac_standard*100}%)", color)
    run_parallel_smina(rec_pdbqt, rapid_surv, standard_out, exh_standard, scoring_func, assigned_cpus, center, size, track_dir)
    count = filter_sdf(standard_out, standard_surv, frac_standard)
    if count == 0: return

    # 3. PRECISE
    precise_out = os.path.join(track_dir, "results_PRECISE.sdf")
    precise_surv = os.path.join(track_dir, "survivors_PRECISE.sdf")
    exh_precise = config.get('exh_precise', 32); frac_precise = config.get('frac_precise', 1.0)
    
    log(scoring_func, f"Step 3: PRECISE (Exh={exh_precise}, Keep={frac_precise*100}%)", color)
    run_parallel_smina(rec_pdbqt, standard_surv, precise_out, exh_precise, scoring_func, assigned_cpus, center, size, track_dir)
    count = filter_sdf(precise_out, precise_surv, frac_precise)
    if count == 0: return
    
    # Generate Complexes
    complex_dir = os.path.join(track_dir, "complexes")
    log(scoring_func, "Generating Receptor-Ligand PDB complexes for MD...", color)
    gen_count = generate_complexes(rec_pdb, precise_surv, complex_dir, num_workers=assigned_cpus)
    log(scoring_func, f"Generated {gen_count} PDB complexes in {complex_dir}", color)

    # 4. Analysis
    analysis_json = os.path.join(track_dir, "analysis_results.json")
    
    run_gbsa = not config.get('no_mmgbsa', False)
    run_plip = not config.get('no_plip', False)
    track_data = analyze_track_results(precise_surv, rec_pdb, track_dir, color, scoring_func, run_gbsa=run_gbsa, run_plip=run_plip)
    
    with open(analysis_json, 'w') as f:
        json.dump(track_data, f, indent=4)
            
    return_dict[scoring_func] = track_data
    log(scoring_func, "Track Finished Successfully.", color)