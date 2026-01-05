import os
import multiprocessing
import argparse
import csv
from pipeline_worker_functions import run_screening_track, ensure_directory
from ligand_prep import run_ligand_prep

# FORMAT GBSA
def format_gbsa(val):
    """Helper to safely format GBSA values."""
    if isinstance(val, (int, float)):
        return f"{val:.2f}"
    return str(val)

# GENERATE REPORTS 
def generate_reports(results_vina, results_vinardo, output_base, scoring_mode):
    """
    Generates two reports:
    1. A Text Table (.txt) for quick reading.
    2. A CSV File (.csv) for Excel/Data Analysis.
    """
    all_ligands = sorted(set(results_vina.keys()) | set(results_vinardo.keys()))
    
    # Define Filenames
    txt_file = os.path.join(output_base, "FINAL_SUMMARY.txt")
    csv_file = os.path.join(output_base, "FINAL_SUMMARY.csv")
    
    # GENERATE CSV
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # CSV Header
        header = [
            'LIGAND', 
            'VINA_SCORE', 'VINARDO_SCORE', 'SCORE_DIFF', 
            'GBSA_VINA', 'GBSA_VINARDO', 'GBSA_DIFF',
            'PLIP_VINA', 'PLIP_VINARDO', 
            'CONSENSUS_HIT' 
        ]
        writer.writerow(header)
        
        for lig in all_ligands:
            v_data = results_vina.get(lig, {'smina':0, 'gbsa':0, 'plip':'-'})
            vo_data = results_vinardo.get(lig, {'smina':0, 'gbsa':0, 'plip':'-'})
            
            # Calculations
            s_v = v_data['smina']
            s_vo = vo_data['smina']
            diff_smina = abs(s_v - s_vo)

            g_v = v_data['gbsa']
            g_vo = vo_data['gbsa']
            
            if isinstance(g_v, (int, float)) and isinstance(g_vo, (int, float)):
                diff_gbsa = abs(g_v - g_vo)
            else:
                diff_gbsa = "N/A"
            
            # Consensus Boolean
            is_consensus = (s_v != 0 and s_vo != 0)
            
            # Write Row
            writer.writerow([
                lig, 
                f"{s_v:.2f}", f"{s_vo:.2f}", f"{diff_smina:.2f}",
                format_gbsa(g_v), format_gbsa(g_vo), format_gbsa(diff_gbsa),
                v_data['plip'], vo_data['plip'],
                str(is_consensus)
            ])
            
    print(f"[REPORT] CSV file generated at: {csv_file}")
    print(f"         (Open this file in Excel to view wide rows without wrapping)")

    # GENERATE TEXT TABLE
    w_lig = len("LIGAND")
    w_plip_v = len("PLIP (VINA)")
    w_plip_vo = len("PLIP (VINARDO)")

    for lig in all_ligands:
        w_lig = max(w_lig, len(lig))
        p_v = str(results_vina.get(lig, {}).get('plip', '-'))
        w_plip_v = max(w_plip_v, len(p_v))
        p_vo = str(results_vinardo.get(lig, {}).get('plip', '-'))
        w_plip_vo = max(w_plip_vo, len(p_vo))
    
    w_lig += 2; w_plip_v += 2; w_plip_vo += 2

    header_str = (
        f"{'LIGAND':<{w_lig}} | "
        f"{'VINA':<7} | {'VINARDO':<7} | {'DIFF(S)':<7} | "
        f"{'GBSA(V)':<8} | {'GBSA(VO)':<8} | {'DIFF(G)':<7} | "
        f"{'PLIP (VINA)':<{w_plip_v}} | {'PLIP (VINARDO)':<{w_plip_vo}}"
    )
    divider = "-" * len(header_str)
    
    consensus_rows = []

    with open(txt_file, 'w') as f:
        f.write(header_str + "\n" + divider + "\n")
        print("\n" + header_str)
        print(divider)
        
        for lig in all_ligands:
            v_data = results_vina.get(lig, {'smina':0, 'gbsa':0, 'plip':'-'})
            vo_data = results_vinardo.get(lig, {'smina':0, 'gbsa':0, 'plip':'-'})
            
            s_v = v_data['smina']
            s_vo = vo_data['smina']
            diff_smina = abs(s_v - s_vo)
            
            g_v = v_data['gbsa']
            g_vo = vo_data['gbsa']
            
            if isinstance(g_v, (int, float)) and isinstance(g_vo, (int, float)):
                d_g = f"{abs(g_v - g_vo):.2f}"
            else:
                d_g = "-"

            row = (
                f"{lig:<{w_lig}} | "
                f"{s_v:<7.2f} | {s_vo:<7.2f} | {diff_smina:<7.2f} | "
                f"{format_gbsa(g_v):<8} | {format_gbsa(g_vo):<8} | {d_g:<7} | "
                f"{str(v_data['plip']):<{w_plip_v}} | {str(vo_data['plip']):<{w_plip_vo}}"
            )
            
            f.write(row + "\n")
            print(row)
            
            if s_v != 0 and s_vo != 0:
                consensus_rows.append(row)

        if consensus_rows:
            cons_header = "\n\n" + "="*40 + " CONSENSUS CANDIDATES " + "="*40
            f.write(cons_header + "\n")
            print(cons_header)
            f.write(header_str + "\n" + divider + "\n")
            print(header_str)
            print(divider)
            for row in consensus_rows:
                f.write(row + "\n")
                print(row)

if __name__ == '__main__':

    multiprocessing.set_start_method('spawn', force=True)
    parser = argparse.ArgumentParser(description="Parallel Virtual Screening Pipeline")
    
    # Required Files
    parser.add_argument('--receptor_pdb', required=True, help='Path to receptor PDB file')
    parser.add_argument('--ligand', required=True, help='Path to ligand SDF file')
    parser.add_argument('--output_dir', default='pipeline_results', help='Directory to save results\n')

    # Box Coordinates
    parser.add_argument('--center_x', type=float, required=True, help='Center X')
    parser.add_argument('--center_y', type=float, required=True, help='Center Y')
    parser.add_argument('--center_z', type=float, required=True, help='Center Z\n')
    parser.add_argument('--size_x', type=float, required=True, help='Size X')
    parser.add_argument('--size_y', type=float, required=True, help='Size Y')
    parser.add_argument('--size_z', type=float, required=True, help='Size Z\n')

    parser.add_argument('--cpu_count', type=int, default=0, help='Total CPUs to allocate (Selecting 0 causes the tool to utilises max. no. of CPUs available)')
    parser.add_argument('--prep_cpus', type=int, default=0, help='CPUs for ligand preparation (0 = Auto/75%%)')

    # Simulation Params
    parser.add_argument('--exh_rapid', type=int, default=2, help='Exhaustiveness RAPID (def: 1)')
    parser.add_argument('--exh_standard', type=int, default=8, help='Exhaustiveness STANDARD (def: 8)')
    parser.add_argument('--exh_precise', type=int, default=32, help='Exhaustiveness PRECISE (def: 32)\n')
    
    parser.add_argument('--frac_rapid', type=float, default=0.5, help='Fraction RAPID (def: 0.5)')
    parser.add_argument('--frac_standard', type=float, default=0.3, help='Fraction STANDARD (def: 0.3)')
    parser.add_argument('--frac_precise', type=float, default=1.0, help='Fraction PRECISE (def: 1.0)\n')
    
    parser.add_argument('--scoring', choices=['vina', 'vinardo', 'both'], default='both', help='Scoring function (def: both)')
    parser.add_argument('--no_mmgbsa', action='store_true', help='Skip MM-GBSA')
    parser.add_argument('--no_plip', action='store_true', help='Skip PLIP Analysis')
    parser.add_argument('--lipinski', action='store_true', help="Filter input ligands using Lipinski's Rule of 5")
    
    args = parser.parse_args()

    max_sys_cpus = multiprocessing.cpu_count()
    user_cpus = args.cpu_count

    if user_cpus > max_sys_cpus or user_cpus < 0:
        print(f"[CONFIGURATION] Requested CPUs ({user_cpus}) invalid or greater than the number of CPUs available ({max_sys_cpus}). Using maximum number of CPUs available.")
        total_cpus = max_sys_cpus
    else:
        total_cpus = user_cpus

    scoring_mode=args.scoring
    if scoring_mode == 'both':
        cpus_per_track = max(1, int(total_cpus//2))
        print(f"[CONFIGURATION] Scoring Mode: BOTH. Allocating {cpus_per_track} CPUs per track (Total: {cpus_per_track*2})")
    else:
        cpus_per_track = total_cpus
        print (f"[CONFIGURATION] Scoring Mode: {scoring_mode.upper()}. Allocating {cpus_per_track} CPUs to track.")

    # SETUP
    output_base = args.output_dir
    ensure_directory(output_base)
    
    config = vars(args)

    config['cpu_count'] = cpus_per_track

    # LIGAND PREPARATION
    print ("===== LIGAND PREPARATION =====")
    prepared_lig_path = run_ligand_prep(
        input_file=args.ligand,
        output_dir=args.output_dir,
        apply_lipinski=args.lipinski,
        cpu_count=args.prep_cpus 
    )
    
    if not prepared_lig_path:
        print("[CRITICAL ERROR] Ligand preparation failed. Exiting....")
        exit(1)

    config['ligand'] = prepared_lig_path
    
    manager = multiprocessing.Manager()
    final_results = manager.dict()
    processes = []

    scoring_mode = args.scoring
    modes_to_run = ['vina', 'vinardo'] if scoring_mode == 'both' else [scoring_mode]
        
    print(f"===== LAUNCHING PIPELINE: {args.scoring.upper()} =====")
    print(f"Receptor: {args.receptor_pdb}")
    print(f"Ligand:   {args.ligand}")
    
    for mode in modes_to_run:
        p = multiprocessing.Process(
            target=run_screening_track, 
            args=(mode, output_base, final_results, config)
        )
        processes.append(p)
        p.start()
        
    for p in processes:
        p.join()
        
    print("\n===== ALL TRACKS COMPLETED. GENERATING REPORTS =====")
    
    # Handle report generation logic based on mode
    if scoring_mode == 'both':
        generate_reports(final_results['vina'], final_results['vinardo'], output_base, 'both')
    elif scoring_mode == 'vina':
        generate_reports(final_results['vina'], {}, output_base, 'vina')
    else:
        generate_reports({}, final_results['vinardo'], output_base, 'vinardo')
        
    print(f"[Done] Pipeline Finished.")