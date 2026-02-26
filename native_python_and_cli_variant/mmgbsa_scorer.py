import sys
import openmm as mm
from openmm import app, unit, Platform, LangevinIntegrator
from openmm.app import Modeller, PDBFile, ForceField, Simulation
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit.topology import Molecule
import time
import os
from pdbfixer import PDBFixer

class GBSAScorer:
    """
    Helper class to perform MM_GBSA rescoring on docked poses.
    Fixed to use NoCutoff for accurate implicit solvent calculations.
    """

    def __init__(self, receptor_pdb_path, forcefield_xmls=None, temp=300.0, friction=1.0, timestep_fs = 2.0, gpu_device_index=0):
        print ("\n===== Performing MM-GBSA =====\n")
        print (f"\n[GBSA] Initializing Force Fields for Receptor.\n")

        if forcefield_xmls is None:
            self.ff_xmls = ['amber14-all.xml', 'implicit/gbn2.xml']
        else:
            self.ff_xmls = forcefield_xmls

        self.temp = temp
        self.friction = friction
        self.timestep_fs = timestep_fs

        self.base_forcefield = ForceField(*self.ff_xmls)

        print(f"[GBSA] Prepping Receptor: {receptor_pdb_path}")
        
        fixer = PDBFixer(filename=receptor_pdb_path)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        modeller = Modeller(fixer.topology, fixer.positions)
        atoms_to_delete = [atom for atom in modeller.topology.atoms() if atom.element.symbol == 'H']
        modeller.delete(atoms_to_delete)
        
        print(f"[GBSA] Adding Hydrogens using Force Field...")
        modeller.addHydrogens(self.base_forcefield, pH=7.0)
        
        base_name = os.path.splitext(receptor_pdb_path)[0]
        fixed_path = f"{base_name}_fixed.pdb"
        
        with open(fixed_path, 'w') as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f)
            
        print(f"[GBSA] Fixed receptor saved to: {fixed_path}")
        self.receptor_pdb = PDBFile(fixed_path)
        self.platform, self.platform_props = self._detect_and_validate_platform(gpu_device_index)

    def _detect_and_validate_platform(self, device_index=0):
        priorities = ['CUDA', 'Metal', 'OpenCL', 'CPU']
        num = Platform.getNumPlatforms()
        available_names = set([Platform.getPlatform(i).getName() for i in range(num)])
        
        for name in priorities:
            if name in available_names:
                try:
                    platform = Platform.getPlatformByName(name)
                    props = {}
                    if name in ['CUDA', 'OpenCL']:
                        props = {'Precision': 'mixed', 'DeviceIndex': str(device_index)}
                    
                    system_dummy = mm.System()
                    system_dummy.addParticle(1.0 * unit.dalton)
                    integrator_dummy = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
                    context_dummy = mm.Context(system_dummy, integrator_dummy, platform, props)
                    del context_dummy, integrator_dummy, system_dummy
                    
                    print(f"[GBSA] Platform '{name}' successfully validated and selected.")
                    return platform, props
                except Exception as e:
                    print(f"[GBSA] Warning: Platform '{name}' failed validation: {e}")
                    continue
        return Platform.getPlatformByName('Reference'), {}

    def calculate_single_dG(self, ligand_mol_openff):
        t_start = time.time()
        print("\n > [1/2] Parameterizing Ligand....", end="", flush=True)
        t1 = time.time()

        forcefield = ForceField(*self.ff_xmls)
        smirnoff = SMIRNOFFTemplateGenerator(molecules=ligand_mol_openff, forcefield='openff-2.1.0')
        forcefield.registerTemplateGenerator(smirnoff.generator)

        print (f"Done ({time.time()-t1:.1f}s).")

        print (" > [2/2] Minimizing Complex....", end="", flush=True)
        t2 = time.time()
        
        modeller = Modeller(self.receptor_pdb.topology, self.receptor_pdb.positions)
        num_rec_atoms = modeller.topology.getNumAtoms()

        lig_top = ligand_mol_openff.to_topology().to_openmm()
        lig_pos = ligand_mol_openff.conformers[0].to_openmm()
        modeller.add(lig_top, lig_pos)

        system = forcefield.createSystem(modeller.topology, 
                                         nonbondedMethod=app.NoCutoff, 
                                         constraints=app.HBonds, 
                                         rigidWater=True)
        
        integrator = LangevinIntegrator(self.temp * unit.kelvin, self.friction / unit.picosecond, self.timestep_fs * unit.femtoseconds)
        simulation = Simulation(modeller.topology, system, integrator, self.platform, self.platform_props)
        simulation.context.setPositions(modeller.positions)
        
        simulation.minimizeEnergy()
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        E_complex = state.getPotentialEnergy()
        all_positions = state.getPositions() 
        
        rec_positions = all_positions[:num_rec_atoms]
        lig_positions = all_positions[num_rec_atoms:]
        print (f"Done ({time.time()-t2:.1f}s).")

        # COMPONENT ENERGIES
        # Receptor
        res_system = self.base_forcefield.createSystem(self.receptor_pdb.topology,
                                                  nonbondedMethod=app.NoCutoff,
                                                  constraints=app.HBonds)
        res_int = LangevinIntegrator(self.temp * unit.kelvin, self.friction / unit.picosecond, self.timestep_fs * unit.femtoseconds)
        res_sim = Simulation(self.receptor_pdb.topology, res_system, res_int, self.platform, self.platform_props)
        res_sim.context.setPositions(rec_positions) 
        E_receptor = res_sim.context.getState(getEnergy=True).getPotentialEnergy()

        # Ligand
        lig_system = forcefield.createSystem(lig_top,
                                             nonbondedMethod=app.NoCutoff,
                                             constraints=app.HBonds)
        lig_int = LangevinIntegrator(self.temp * unit.kelvin, self.friction / unit.picosecond, self.timestep_fs * unit.femtoseconds)
        lig_sim = Simulation(lig_top, lig_system, lig_int, self.platform, self.platform_props)
        lig_sim.context.setPositions(lig_positions)
        E_ligand = lig_sim.context.getState(getEnergy=True).getPotentialEnergy()

        # Convert to kcal/mol
        dG = (E_complex - (E_receptor + E_ligand)).value_in_unit(unit.kilocalories_per_mole)
        
        e_c_kcal = E_complex.value_in_unit(unit.kilocalories_per_mole)
        e_r_kcal = E_receptor.value_in_unit(unit.kilocalories_per_mole)
        e_l_kcal = E_ligand.value_in_unit(unit.kilocalories_per_mole)
        print(f"E_complex: {e_c_kcal:.1f} | E_receptor: {e_r_kcal:.1f} | E_ligand: {e_l_kcal:.1f} | dG: {dG:.1f}")

        total_time = time.time() - t_start
        print (f" > [Total time : {total_time:.1f}s]")

        return dG
    
    def process_sdf_file(self, sdf_path):
        results = []
        print(f"[GBSA] Reading ligands from {sdf_path}....")
        try:
            molecules = Molecule.from_file(sdf_path, allow_undefined_stereo=True)
        except Exception as e:
            print(f"[ERROR] Failed to load SDF: {e}")
            return []
        
        if isinstance(molecules, Molecule): molecules = [molecules]
        print(f"[GBSA] Found {len(molecules)} candidates.")

        for i, mol in enumerate(molecules):
            try:
                name = mol.name if mol.name else f"Ligand_{i+1}"
                smina_score = mol.properties.get('minimizedAffinity', "N/A")
                print(f"Processing {name}....")
                gbsa_score = self.calculate_single_dG(mol)
                results.append({"Name": name, "Binding Affinity Score": smina_score, "GBSA Score": gbsa_score})
            except Exception as e:
                print(f"[ERROR] Skipping {name}: {e}")
                results.append({"Name": name, "Binding Affinity Score": smina_score, "GBSA Score": "ERROR"})
        return results