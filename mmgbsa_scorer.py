import sys
from openmm import app, unit, Platform, LangevinIntegrator
from openmm.app import Modeller, PDBFile, ForceField, Simulation
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit.topology import Molecule
import time

class GBSAScorer:
    """
    Helper class to perform MM_GBSA rescoring on docked poses.
    """

    def __init__(self, receptor_pdb_path):
        """
        Initializing with receptor PDB + ForceField loading for receptor.
        """
        print ("\n===== Performing MM-GBSA =====\n")
        print (f"\n[GBSA] Initializing Force Fields for Receptor.\n")

        self.receptor_pdb = PDBFile(receptor_pdb_path)
        self.base_forcefield = ForceField ('amber14-all.xml', 'implicit/obc2.xml')

        self.platform, self.platform_props = self._detect_and_validate_platform()

    def _detect_and_validate_platform(self):
        """
        Iterates through preferred platforms (CUDA > Metal > OpenCL > CPU).
        Validates them by attempting to create a dummy simulation context.
        If a platform exists but crashes (e.g. broken OpenCL driver), it skips it.
        """
        # Preference Order
        priorities = ['CUDA', 'Metal', 'OpenCL', 'CPU']
        
        num = Platform.getNumPlatforms()
        available_names = set([Platform.getPlatform(i).getName() for i in range(num)])
        
        for name in priorities:
            if name in available_names:
                try:
                    # 1. Configuration
                    platform = Platform.getPlatformByName(name)
                    props = {}
                    if name in ['CUDA', 'OpenCL']:
                        props = {'Precision': 'mixed'}
                    
                    # 2. VALIDATION: Try to launch a dummy simulation
                    system_dummy = app.System()
                    system_dummy.addParticle(1.0 * unit.dalton)
                    integrator_dummy = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
                    
                    context_dummy = app.Context(system_dummy, integrator_dummy, platform, props)
                    
                    del context_dummy, integrator_dummy, system_dummy
                    
                    print(f"[GBSA] Platform '{name}' successfully validated and selected.")
                    return platform, props

                except Exception as e:
                    print(f"[GBSA] Warning: Platform '{name}' is listed but failed validation. Skipping.\n       Error: {e}")
                    continue
        
        print("[GBSA] Critical Warning: No functional platform found. Defaulting to Reference (Slowest).")
        return Platform.getPlatformByName('Reference'), {}

    def calculate_single_dG(self, ligand_mol_openff):
        t_start = time.time()

        print(" > [1/4] Parameterizing Ligand (OpenFF Sage)....", end="", flush=True)
        t1 = time.time()

        forcefield = ForceField ('amber14-all.xml', 'implicit/obc2.xml')

        smirnoff = SMIRNOFFTemplateGenerator(molecules=ligand_mol_openff, forcefield='openff-2.1.0')
        forcefield.registerTemplateGenerator(smirnoff.generator)
        print (f"Done ({time.time()-t1:.1f}s).")

        print (" > [2/4] Building System....", end="", flush=True)
        t2 = time.time()
        
        # Create Modeller
        modeller = Modeller(self.receptor_pdb.topology, self.receptor_pdb.positions)
        
        # Add Ligand
        lig_top = ligand_mol_openff.to_topology().to_openmm()
        lig_pos = ligand_mol_openff.conformers[0].to_openmm()
        modeller.add(lig_top, lig_pos)

        # Create System
        system = forcefield.createSystem(modeller.topology, 
                                         nonbondedMethod=app.CutoffNonPeriodic, 
                                         nonbondedCutoff=2.0*unit.nanometers,
                                         constraints=app.HBonds, 
                                         rigidWater=True)
        
        integrator = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)

        # Use the pre-validated platform
        simulation = Simulation(modeller.topology, system, integrator, self.platform, self.platform_props)
        simulation.context.setPositions(modeller.positions)
        print (f"Done ({time.time()-t2:.1f}s).")

        print (" > [3/4] Minimizing Energy....", end="", flush=True)
        t3 = time.time()
        simulation.minimizeEnergy()
        print (f"Done ({time.time()-t3:.1f}s).")

        print (" > [4/4] Calculating Component Energies....", end="", flush=True)
        
        # Calculate Energies
        # 1. Complex Energy
        state = simulation.context.getState(getEnergy=True)
        E_complex = state.getPotentialEnergy()

        # 2. Receptor Energy (Remove ligand atoms)
        res_system = self.base_forcefield.createSystem(self.receptor_pdb.topology,
                                                  nonbondedMethod=app.CutoffNonPeriodic,
                                                  nonbondedCutoff=2.0*unit.nanometers,
                                                  constraints=app.HBonds)
        res_int = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
        res_sim = Simulation(self.receptor_pdb.topology, res_system, res_int, self.platform, self.platform_props)
        res_sim.context.setPositions(self.receptor_pdb.positions)
        E_receptor = res_sim.context.getState(getEnergy=True).getPotentialEnergy()

        # 3. Ligand Energy
        lig_system = forcefield.createSystem(lig_top,
                                             nonbondedMethod=app.CutoffNonPeriodic,
                                             nonbondedCutoff=2.0*unit.nanometers,
                                             constraints=app.HBonds)
        lig_int = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
        lig_sim = Simulation(lig_top, lig_system, lig_int, self.platform, self.platform_props)
        lig_sim.context.setPositions(lig_pos)
        E_ligand = lig_sim.context.getState(getEnergy=True).getPotentialEnergy()

        # FINAL DELTA G
        dG = E_complex - (E_receptor + E_ligand)
        total_time = time.time() - t_start
        print (f" > [Total time : {total_time:.1f}s]")

        return dG.value_in_unit(unit.kilocalories_per_mole)
    
    
    def process_sdf_file(self, sdf_path):
        """
        Reads an SDF, calculates MM-GBSA for all the ligands in the file, and returns results.
        """
        results = []
        print(f"[GBSA] Reading ligands from {sdf_path}....")

        try:
            molecules = Molecule.from_file(sdf_path, allow_undefined_stereo=True)
        except Exception as e:
            print(f"[ERROR] Failed to load SDF: {e}")
            return []
        
        # Handle Single Molecule Case
        if isinstance(molecules, Molecule):
            molecules = [molecules]
        
        print(f"[GBSA] Found {len(molecules)} candidates. Starting MM-GBSA Calculation....")

        for i, mol in enumerate(molecules):
            try:
                name = mol.name if mol.name else f"Ligand_{i+1}"
                smina_score = mol.properties.get('minimizedAffinity', "N/A")
                print(f"Processing {name}....")

                gbsa_score = self.calculate_single_dG(mol)

                results.append({
                    "Name": name,
                    "Binding Affinity Score": smina_score,
                    "GBSA Score": gbsa_score
                })

            except Exception as e:
                print(f"[ERROR] Skipping {name} due to error: {e}")
                results.append({
                    "Name": name,
                    "Binding Affinity Score": smina_score,
                    "GBSA Score": "ERROR"
                })

        return results