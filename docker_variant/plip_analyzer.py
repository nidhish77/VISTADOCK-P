import os
import subprocess
import xml.etree.ElementTree as ET
from openmm.app import PDBFile, Modeller
from openff.toolkit.topology import Molecule
from collections import defaultdict

class PLIPAnalyzer:
    def __init__(self, receptor_pdb_path):
        self.receptor_path = receptor_pdb_path
        self.receptor_pdb = PDBFile(receptor_pdb_path)

    def analyze_interactions(self, ligand_mol, output_dir="plip_reports"):
        """
        Creates a complex, runs PLIP and returns summary of interactions.
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        lig_name = ligand_mol.name if ligand_mol.name else "Ligand"
        lig_name = "".join([c for c in lig_name if c.isalnum() or c in ('_','-')])
        complex_filename = os.path.join(output_dir, f"complex_{lig_name}.pdb")

        try:
            modeller = Modeller(self.receptor_pdb.topology, self.receptor_pdb.positions)
            lig_top = ligand_mol.to_topology().to_openmm()
            lig_pos = ligand_mol.conformers[0].to_openmm()

            modeller.add(lig_top, lig_pos)

            with open (complex_filename, 'w') as f:
                PDBFile.writeFile(modeller.topology, modeller.positions, f)
        
        except Exception as e:
            return f"[Error building complex: {e}]"
        
        cmd = [
            "plip",
            "-f", complex_filename,     #Input file
            "-x",                       #Output XML
            "-o", output_dir,           #Output directory
            "--quiet"                   #Reduced terminal noise
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"\n[PLIP ERROR] Command failed for {lig_name}:")
            print(f"    STDOUT: {e.stdout}")
            print(f"    STDERR: {e.stderr}")
            return "PLIP execution failed."
        except FileNotFoundError:
            print("\n[PLIP ERROR] 'plip' command not found. Please install it: 'conda install -c conda-forge plip'")
            return "PLIP not found."
        
        expected_xml = os.path.join(output_dir, "report.xml")
        if not os.path.exists(expected_xml):
            # Fallback: check for files named with the complex prefix or just the newest XML
            import glob
            xmls = glob.glob(os.path.join(output_dir, "*.xml"))
            if xmls:
                expected_xml = max(xmls, key=os.path.getctime)
            else:
                return "No PLIP XML found."
        
        return self.parse_plip_xml(expected_xml)
    
    def parse_plip_xml(self, xml_file):
        """
        Extracts interactions and residues from XML file.
        Returns a grouped string: 'Type: Res1, Res2 | Type: Res3'
        """
        interactions = defaultdict(set)
        
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()

            for site in root.findall(".//bindingsite"):
                # Hydrogen Bonds
                for hbond in site.findall(".//hydrogen_bond"):
                    res = hbond.find("resnr").text
                    restype = hbond.find("restype").text
                    interactions['HBond'].add(f"{restype}{res}")

                # Hydrophobic Interactions
                for hydro in site.findall(".//hydrophobic_interaction"):
                    res = hydro.find("resnr").text
                    restype = hydro.find("restype").text
                    interactions['Hydrophobic'].add(f"{restype}{res}")

                # Salt Bridges
                for salt in site.findall(".//salt_bridge"):
                    res = salt.find("resnr").text
                    restype = salt.find("restype").text
                    interactions['Salt Bridge'].add(f"{restype}{res}")

                # Pi-Stacking
                for stack in site.findall(".//pi_stacking"):
                    res = stack.find("resnr").text
                    restype = stack.find("restype").text
                    interactions['PiStack'].add(f"{restype}{res}")

                # Halogen Bonds
                for halogen in site.findall(".//halogen_bond"):
                    res = halogen.find("resnr").text
                    restype = halogen.find("restype").text
                    interactions['Halogen'].add(f"{restype}{res}")

            if not interactions:
                return "No Specific Interactions"

            # Output Formatting: "Type: Res1, Res2 | Type: Res3"
            output_parts = []
            for itype in sorted(interactions.keys()):
                residues = sorted(list(interactions[itype]))
                res_str = ", ".join(residues)
                output_parts.append(f"{itype}: {res_str}")
            
            return " | ".join(output_parts)
        
        except Exception as e:
            return f"XML parse error: {e}"

    def clean_plip_files(output_dir):
        import shutil
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)