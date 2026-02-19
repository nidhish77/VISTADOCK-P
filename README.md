# VISTADOCK-P

**VISTADOCK-P** is a high-throughput virtual screening pipeline designed to prioritize high-affinity ligands through a *funnel-like* protocol. It integrates **Consensus Docking** (Vina & Vinardo scoring algorithms via SMINA), **Binding Free Energy Re-scoring** (MM-GBSA via OpenMM), **Convolutional Neural Network Re-scoring** (via GNINA) and **Residue Interaction Analysis** (via PLIP) into an automated loop for effectively screening thousands of ligands.

This repository contains **two** variants of the pipeline:
1. **Docker Variant**: Runs in an isolated Docker container (*recommended for stability*)
2. **Native Python Variant**: Runs directly on your host machine (*recommended for modularity and advanced users*)

------------------------------------------------------------------------
## Option 1: DOCKER VARIANT
Located in the *`/docker_variant`* folder.

This version is pre-packaged with all the dependencies required to run the pipeline.

### 1. Pre-requisites:
- Docker Desktop
- NVIDIA Drivers (for CNN Re-scoring and faster MM-GBSA calculations)

### 2. Setup:
Navigate to the *`/docker_variant`* directory and build the Docker image.

```bash
cd docker_variant
sudo docker build --platform linux/amd64 -t vistadock_p .
```

The pipeline has a GUI that can be used to run it. The GUI Launcher is packaged within the Docker image, so if required, the `gui_launcher.py` script can be extracted from the image using the following command.

```bash
sudo docker run --rm vistadock_p cat /app/gui_launcher.py > gui_launcher.py
```

### 3. Running the Pipeline:
Open the GUI Launcher using the following command.
```bash
python3 gui_launcher.py
```
Read the Disclaimers within the GUI, input all your values in the provided fields, then click on the `GENERATE COMMAND` button. Once the command has been generated in the text box below the button, click on the `RUN PIPELINE` button, and the pipeline will automatically run on your terminal. 

Alternatively you can copy the command that is generated, and run it directly on your command line.

------------------------------------------------------------------------
## Option 2: NATIVE VARIANT
Located in the *`/native_python_and_cli_variant`* folder.

This version directly on your OS via a Conda environment.

### 1. Pre-requisites:
- Anaconda / Miniconda / Mamba
- `gnina` binary installed in your system PATH.

### 2. Setup:
Navigate to the *`/native_python_and_cli_variant`* directory and create the environment.

```bash
cd native_python_and_cli_variant
conda env create -f environment.yml
conda activate vistadock_p_env
```

<<<<<<< HEAD
Download the GNINA binary as follows:
```bash
wget https://github.com/gnina/gnina/releases/download/v1.0.3/gnina
chmod +x gnina
mv gnina $CONDA_PREFIX/bin/  #To ensure GNINA only affects your conda environment
```

=======
>>>>>>> fda99793d70333faa907462fc031f713b1755617
### 3. Running the Pipeline:
You can run the pipeline using the GUI Launcher script present in the *`/native_python_and_cli_variant`* folder.

```bash
python3 gui_launcher.py
```

Read the Disclaimers within the GUI, input all your values in the provided fields, then click on the `GENERATE COMMAND` button. Once the command has been generated in the text box below the button, click on the `RUN PIPELINE` button, and the pipeline will automatically run on your terminal. 

Alternatively you can copy the command that is generated, and run it directly on your command line.
```bash
python3 -u vistadock.py --receptor_pdb receptor.pdb --ligand ligands.sdf ....
```

------------------------------------------------------------------------

## How it works (The 'Survivor-based' *Funnel* Protocol)

VISTADOCK-P processes ligands in three stages of increasing exhaustiveness: 

1. **RAPID Stage**: Fast sweep through all the available ligands at low exhaustiveness *(default: 2)* and only specified fraction of ligands *(default: 0.5)* survive to next stage.

2. **BALANCED Stage**: Higher Precision re-docking of survivors from the RAPID Stage at a higher exhaustiveness *(default: 8)*. Only specified fraction of ligands *(default: 0.3)* survive to next stage.

3. **ULTRA Stage**: Intensive sampling of survivors from the BALANCED Stage at a very high exhaustiveness *(default: 32)*. Specified fraction of survivors from this stage *(default: 1.0)* are carried on to the CNN Re-scoring / Binding Free Energy Re-scoring stage.

4. **CNN Re-scoring** (if selected and GPU available): Survivors from the ULTRA Stage are re-scored using a Convolutional Neural Network via **GNINA**.

5. **Binding Free Energy Re-scoring** (if selected): Receptor-Ligand complexes from the ULTRA Stage are re-scored using their $\Delta$G values via **MM-GBSA** (using OpenMM). 
	(**NOTE**: The *Energy Minimisation* step of MM-GBSA can be performed using GPU acceleration to increase speed significantly, but if an NVIDIA Device is not available, MM-GBSA can still be performed using only the CPU, albeit at an extremely slow pace.)

6. **Residue Interaction Analysis** (if selected): Receptor-Ligand complexes from the ULTRA Stage undergo residue interaction analysis via **PLIP**.

------------------------------------------------------------------------

## Output Files

- `FINAL_SUMMARY.CSV`: Comprehensive spreadsheet containing scores from each step, as well as consensus booleans (whether a particular ligand was scored by both *Vina* and *Vinardo*).
- `FINAL_SUMMARY.TXT`: Summary table in text format.
- `FINAL_SUMMARY.HTML`: Summary table in HTML format, with links to download the Receptor-Ligand Complex for a specific top-ranked Ligand.
- *`*_track/`*: Subfolders containing the Receptor-Ligand PDB complexes and SDF Files (results & survivors) for every stage.








