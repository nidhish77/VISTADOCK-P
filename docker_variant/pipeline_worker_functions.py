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
import sys
import time

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
    print(f"{color}[{tag.upper()}] {message}{Colors.ENDC}", flush=True)

def ensure_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)

# GPU Detection 
def is_gpu_available():
    if shutil.which("nvidia-smi") is None:
        return False
    try:
        subprocess.check_call("nvidia-smi", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except:
        return False

# Count Molecules 
def count_molecules(sdf_path):
    if not os.path.exists(sdf_path): return 0
    try:
        with open(sdf_path, 'r', errors='ignore') as f: 
            return f.read().count("$$$$")
    except:
        return 0

# Split SDF 
def split_sdf_file(input_sdf, chunks_dir, num_chunks):
    if not os.path.exists(chunks_dir):
        os.makedirs(chunks_dir)
    
    with open(input_sdf, 'r') as f:
        content = f.read()

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
            for mol_block in subset: f.write(mol_block)
        chunk_files.append(chunk_name)
    return chunk_files

def merge_sdfs(chunk_files, output_file):
    with open(output_file, 'w') as outfile:
        for fname in chunk_files:
            if os.path.exists(fname):
                with open(fname, 'r') as infile:
                    outfile.write(infile.read())

def parse_scores_text(filepath):
    """
    Parses a file (SDF or Mol2) and returns a dictionary:
    { 'MoleculeName': { 'TagName': 'TagValue', ... } }
    """
    scores = {}
    ext = os.path.splitext(filepath)[1].lower()
    
    with open(filepath, 'r') as f:
        content = f.read()

    if ext == '.sdf':
        # SDF Parsing
        blocks = content.split("$$$$")
        for block in blocks:
            if not block.strip(): continue
            lines = block.strip().split('\n')
            name = lines[0].strip() # First line is name in SDF
            
            mol_props = {}
            target_tags = ['CNNscore', 'CNNaffinity', 'CNN_VS', 'CNNaffinity_variance', 'minimizedAffinity']
            
            for tag in target_tags:
                tag_marker = f"> <{tag}>"
                if tag_marker in block:
                    try:
                        start = block.find(tag_marker)
                        sub = block[start:]
                        val_line = sub.split('\n')[1].strip()
                        mol_props[tag] = val_line
                    except:
                        pass
            
            if name and mol_props:
                scores[name] = mol_props

    elif ext == '.mol2':
        delimiter = "@<TRIPOS>MOLECULE"
        blocks = content.split(delimiter)
        for block in blocks:
            if not block.strip(): continue
            lines = block.strip().split('\n')

            name = lines[0].strip()
            
            mol_props = {}

            for line in lines:
                if "CNNscore" in line or "CNNaffinity" in line or "minimizedAffinity" in line:
                    parts = line.replace('#', '').strip().split()
                    if len(parts) >= 2:
                        for t in ['CNNscore', 'CNNaffinity', 'CNN_VS', 'minimizedAffinity']:
                            if t in parts[0]:
                                mol_props[t] = parts[-1]
            
            if name and mol_props:
                scores[name] = mol_props

    return scores

def merge_cnn_scores(clean_file, scored_file, output_file):
    """
    Grafts scores from scored_file onto clean_file using pure text manipulation.
    Works for SDF and Mol2.
    """
    if not os.path.exists(clean_file) or not os.path.exists(scored_file):
        return False

    # Extract Scores
    score_map = parse_scores_text(scored_file)
    if not score_map:
        print("[WARNING] No scores found in scored file to merge.")
        return False

    ext = os.path.splitext(clean_file)[1].lower()
    
    with open(clean_file, 'r') as f:
        clean_content = f.read()

    with open(output_file, 'w') as out:
        
        if ext == '.sdf':
            blocks = clean_content.split("$$$$")
            for block in blocks:
                if not block.strip(): continue
                
                # Identify Name
                lines = block.strip().split('\n')
                name = lines[0].strip()
                
                new_block = block.strip()
                
                if name in score_map:
                    # Append tags for SDF 
                    props = score_map[name]
                    for key, val in props.items():
                        # Remove existing tag if present to avoid duplicates
                        if f"> <{key}>" not in new_block:
                            new_block += f"\n\n> <{key}>\n{val}"
                
                out.write(new_block + "\n\n$$$$\n")

        elif ext == '.mol2':
            delimiter = "@<TRIPOS>MOLECULE"
            blocks = clean_content.split(delimiter)
            
            for i, block in enumerate(blocks):
                if not block.strip(): continue
                
                lines = block.strip().split('\n')
                name = lines[0].strip()
                
                new_block = block
                
                if name in score_map:
                    props = score_map[name]
                    info_str = "\n@<TRIPOS>COMMENT\n"
                    for key, val in props.items():
                        info_str += f"{key}: {val}\n"
                    
                    new_block += info_str
                
                out.write(delimiter + new_block)

    return True

# Pure CPU Docking Command 
def run_docking_cmd(receptor, ligand, output, exhaustiveness, scoring_func, cpu_count, center, size):
    """
    Standard docking command. Strictly CPU.
    """
    gnina_exe = "gnina"
    
    cmd = [
        'smina', '-r', receptor, '-l', ligand, '-o', output,
        '--exhaustiveness', str(exhaustiveness), '--scoring', scoring_func,
        '--cpu', str(cpu_count),
        '--center_x', str(center[0]), '--center_y', str(center[1]), '--center_z', str(center[2]),
        '--size_x', str(size[0]), '--size_y', str(size[1]), '--size_z', str(size[2]),
        '--quiet' 
    ]
    
    # Hide GPU to prevent accidental initialization
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = ""

    try:
        subprocess.run(cmd, check=True, env=env)
        return True
    except subprocess.CalledProcessError:
        return False

# DOCKING: Parallel Engine (CPU Only) 
def run_parallel_docking(receptor, input_sdf, output_sdf, exhaustiveness, scoring_func, total_cpus, center, size, work_dir):
    
    num_mols = count_molecules(input_sdf)
    
    # Threading Mode
    if num_mols < 4 or total_cpus <= 1:
        log(scoring_func, f"   > Mode: INTERNAL THREADING ({num_mols} ligands, {total_cpus} CPUs)", Colors.BLUE)
        run_docking_cmd(receptor, input_sdf, output_sdf, exhaustiveness, scoring_func, total_cpus, center, size)
        return

    # Parallel Process Mode
    log(scoring_func, f"   > Mode: PARALLEL SPLIT ({num_mols} ligands over {total_cpus} workers)", Colors.BLUE)
    
    chunks_dir = os.path.join(work_dir, "chunks_in")
    out_dir = os.path.join(work_dir, "chunks_out")
    if not os.path.exists(out_dir): os.makedirs(out_dir)

    input_chunks = split_sdf_file(input_sdf, chunks_dir, total_cpus)
    if not input_chunks: return 

    processes = []
    output_chunks = []
    
    # Env to hide GPU
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = ""
    
    docking_exe = "smina" 

    for i, chunk_in in enumerate(input_chunks):
        chunk_out = os.path.join(out_dir, f"out_{i}.sdf")
        output_chunks.append(chunk_out)
        
        cmd = [
            docking_exe, '-r', receptor, '-l', chunk_in, '-o', chunk_out,
            '--exhaustiveness', str(exhaustiveness), '--scoring', scoring_func,
            '--cpu', '1', 
            '--center_x', str(center[0]), '--center_y', str(center[1]), '--center_z', str(center[2]),
            '--size_x', str(size[0]), '--size_y', str(size[1]), '--size_z', str(size[2]),
            '--quiet'
        ]

        # Silence the individual workers 
        p = subprocess.Popen(cmd, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        processes.append(p)

    # Main Process monitors the workers and draws a clean bar
    total_procs = len(processes)
    while True:
        # Count how many processes have a return code 
        completed = sum(1 for p in processes if p.poll() is not None)
        
        # Calculate percentage
        percent = int((completed / total_procs) * 100)
        
        # Draw Progress Bar
        bar_len = 30
        filled_len = int(bar_len * completed // total_procs)
        bar = '=' * filled_len + '-' * (bar_len - filled_len)
        
        # Printed with \r to overwrite the line
        sys.stdout.write(f"\r    [{bar}] {percent}% ({completed}/{total_procs} chunks)")
        sys.stdout.flush()
        
        if completed == total_procs:
            break
            
        time.sleep(2.0) # Update 1 time every 2 seconds

    print()

    merge_sdfs(output_chunks, output_sdf)
    shutil.rmtree(chunks_dir, ignore_errors=True)
    shutil.rmtree(out_dir, ignore_errors=True)

def perform_cnn_rescoring(receptor, input_sdf, output_sdf, gpu_id=0):
    """
    Runs a SINGLE Gnina process to rescore all ligands in the SDF using the GPU.
    This prevents Out Of Memory Error because only one model is loaded.
    """
    if not is_gpu_available():
        log("CNN", "GPU not found. Skipping CNN Rescoring.", Colors.WARNING)
        shutil.copy(input_sdf, output_sdf) # Just copy input to output
        return

    gnina_exe = "gnina"
    log("CNN", "Starting Batch CNN Rescoring on GPU...", Colors.GREEN)
    
    # --score_only calculates the score for the existing pose
    # --cnn_scoring adds the CNN tags
    cmd = [
        gnina_exe, '-r', receptor, '-l', input_sdf, '-o', output_sdf,
        '--score_only', '--cnn_scoring', 'rescore', '--device', str(gpu_id), '--quiet'
    ]

    try:
        subprocess.run(cmd, check=True)
        log("CNN", "CNN Rescoring Complete.", Colors.GREEN)
        return True
    except subprocess.CalledProcessError as e:
        log("CNN", f"CNN Rescoring Failed: {e}", Colors.FAIL)
        shutil.copy(input_sdf, output_sdf)
        return False


# Complex Generator 
def _generate_complex_worker(data):
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

# Complex Generation 
def generate_complexes(rec_pdb_path, lig_sdf_path, output_dir, num_workers=1):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    try:
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

def analyze_track_results(ultra_survivors_path, rec_pdb_path, output_dir, tag_color, tag_name, run_gbsa=True, run_plip=True, forcefield_xmls=None, gbsa_params=None, gpu_id=0):
    results = {}
    gb_scorer = GBSAScorer(rec_pdb_path, forcefield_xmls=forcefield_xmls, 
                           temp=gbsa_params.get('temp', 300.0), 
                           friction=gbsa_params.get('friction', 1.0), 
                           timestep_fs=gbsa_params.get('timestep', 2.0),
                           gpu_device_index=gpu_id) if run_gbsa else None
                           
    plip_tool = PLIPAnalyzer(rec_pdb_path) if run_plip else None
    plip_out_dir = os.path.join(output_dir, "plip_reports")
    
    try:
        mols = Molecule.from_file(ultra_survivors_path, allow_undefined_stereo=True)
        if isinstance(mols, Molecule):
            mols = [mols]
    except:
        log(tag_name, "Error loading survivors for analysis.", Colors.FAIL)
        return {}

    msg = f"Analyzing {len(mols)} candidates..." if run_gbsa else f"Profiling interactions for {len(mols)} candidates (GBSA Skipped)..."
    log(tag_name, msg, tag_color)

    for i, mol in enumerate(mols):
        name = mol.name if mol.name else f"Ligand_{i+1}"
        smina_raw = mol.properties.get('minimizedAffinity', "N/A")
        if isinstance(smina_raw, str) and smina_raw != 'N/A':
            smina_raw = smina_raw.strip().split('\n')[0]

        try:
            final_smina = float(smina_raw)
        except ValueError:
            final_smina = 0.0
            
        gbsa = "Skipped"
        if run_gbsa and gb_scorer:
            try: gbsa = gb_scorer.calculate_single_dG(mol)
            except Exception as e: 
                print(f"[{tag_name}] GBSA Error for {name}: {e}") 
                gbsa = -999.9 
            
        interactions = "Skipped"
        if run_plip and plip_tool:
            try: interactions = plip_tool.analyze_interactions(mol, output_dir=plip_out_dir)
            except: interactions = "Error"
            
        results[name] = {
            'smina': float(final_smina) if final_smina != "N/A" else 0.0, 
            'gbsa': gbsa, 
            'plip': interactions,
            'properties': mol.properties
        }
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
    assigned_gpu = config.get('gpu_device', 0)
    log(scoring_func, f"Track Initiated using {assigned_cpus} CPUs.", color)

    # Prep Receptor
    rec_pdbqt = os.path.join(track_dir, "receptor_prepared.pdbqt")
    if not os.path.exists(rec_pdbqt):
        log(scoring_func, "Converting Receptor PDB to PDBQT...", color)
        if not convert_pdb_to_pdbqt(rec_pdb, rec_pdbqt):
            log(scoring_func, "Failed to convert Receptor. Aborting.", Colors.FAIL); return

    # RAPID 
    rapid_out = os.path.join(track_dir, "results_RAPID.sdf")
    rapid_surv = os.path.join(track_dir, "survivors_RAPID.sdf")
    exh_rapid = config.get('exh_rapid', 2); frac_rapid = config.get('frac_rapid', 0.5)
    
    log(scoring_func, f"Step 1: RAPID (Exh={exh_rapid}, Keep={frac_rapid*100}%)", color)
    run_parallel_docking(rec_pdbqt, lig_sdf, rapid_out, exh_rapid, scoring_func, assigned_cpus, center, size, track_dir)
    count = filter_sdf(rapid_out, rapid_surv, frac_rapid)
    if count == 0: return

    # BALANCED 
    balanced_out = os.path.join(track_dir, "results_BALANCED.sdf")
    balanced_surv = os.path.join(track_dir, "survivors_BALANCED.sdf")
    exh_balanced = config.get('exh_balanced', 8); frac_balanced = config.get('frac_balanced', 0.3)
    
    log(scoring_func, f"Step 2: BALANCED (Exh={exh_balanced}, Keep={frac_balanced*100}%)", color)
    run_parallel_docking(rec_pdbqt, rapid_surv, balanced_out, exh_balanced, scoring_func, assigned_cpus, center, size, track_dir)
    count = filter_sdf(balanced_out, balanced_surv, frac_balanced)
    if count == 0: return

    # ULTRA
    ultra_out = os.path.join(track_dir, "results_ULTRA.sdf")
    ultra_surv = os.path.join(track_dir, "survivors_ULTRA.sdf")
    exh_ultra = config.get('exh_ultra', 32); frac_ultra = config.get('frac_ultra', 1.0)
    
    log(scoring_func, f"Step 3: ULTRA (Exh={exh_ultra}, Keep={frac_ultra*100}%)", color)
    run_parallel_docking(rec_pdbqt, balanced_surv, ultra_out, exh_ultra, scoring_func, assigned_cpus, center, size, track_dir)
    count = filter_sdf(ultra_out, ultra_surv, frac_ultra)
    if count == 0: return
    
    # CNN RESCORING (GPU - Single Process)

    final_ligands_path = ultra_surv
    
    final_ligands_path = ultra_surv
    
    if config.get('cnn_scoring', False):
        cnn_raw_out = os.path.join(track_dir, "results_CNN_RAW.sdf") 
        cnn_merged_out = os.path.join(track_dir, "results_CNN_MERGED.sdf")
        
        log(scoring_func, "Step 4: CNN Rescoring (GPU Accelerated)...", color)
        
        # Run GNINA (Outputting to RAW file)
        success = perform_cnn_rescoring(rec_pdbqt, ultra_surv, cnn_raw_out, gpu_id=assigned_gpu)
        
        # Merge Scores 
        if success and os.path.exists(cnn_raw_out):
            log(scoring_func, "   > Grafting CNN scores onto clean topology...", color)
            
            merge_success = merge_cnn_scores(ultra_surv, cnn_raw_out, cnn_merged_out)
            
            if merge_success:
                final_ligands_path = cnn_merged_out
            else:
                log(scoring_func, "   > Merge failed. Using original survivors.", Colors.WARNING)
        else:
             log(scoring_func, "   > CNN Rescoring failed. Proceeding with standard scores.", Colors.WARNING)
            

    # Generate Complexes (uses final_ligands_path)
    complex_dir = os.path.join(track_dir, "complexes")
    log(scoring_func, "Generating Receptor-Ligand PDB complexes for MD...", color)
    gen_count = generate_complexes(rec_pdb, final_ligands_path, complex_dir, num_workers=assigned_cpus)
    log(scoring_func, f"Generated {gen_count} PDB complexes in {complex_dir}", color)

    # Analysis
    analysis_json = os.path.join(track_dir, "analysis_results.json")
    
    run_gbsa = not config.get('no_mmgbsa', False)
    run_plip = not config.get('no_plip', False)
    ff_xmls = config.get('ff_xmls', ['amber14-all.xml', 'implicit/gbn2.xml'])
    
    gbsa_params = {
        'temp': config.get('gbsa_temp', 300.0),
        'friction': config.get('gbsa_friction', 1.0),
        'timestep': config.get('gbsa_timestep', 2.0)
    }

    # Pass the RESCORED file to analysis
    track_data = analyze_track_results(final_ligands_path, rec_pdb, track_dir, color, scoring_func, run_gbsa=run_gbsa, run_plip=run_plip, forcefield_xmls=ff_xmls, gbsa_params=gbsa_params, gpu_id=assigned_gpu)
    
    with open(analysis_json, 'w') as f:
        json.dump(track_data, f, indent=4)
            
    return_dict[scoring_func] = track_data
    log(scoring_func, "Track Finished Successfully.", color)