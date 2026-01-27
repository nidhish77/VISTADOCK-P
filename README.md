# VISTADOCK-P

**VISTADOCK-P** is a high-throughput virtual screening pipeline designed to prioritize high-affinity ligands through a *funnel-like* protocol. It integrates **Consensus Docking** (Vina & Vinardo scoring algorithms via SMINA), **Binding Free Energy Re-scoring** (MM-GBSA via OpenMM), **Convolutional Neural Network Re-scoring** (via GNINA) and  **Interaction Topology** (via PLIP) into an automated loop for effectively screening thousands of ligands.

This repository contains **two** variants of the pipeline:
1. **Docker Variant**: Runs in an isolated Docker container (*recommended for stability*)
2. **Native Python Variant**: Runs directly on your host machine (*recommended for modularity and advanced users*)

## Option 1: DOCKER VARIANT
Located in the *`\docker_variant`* folder.

This version is pre-packaged with all the dependencies required to run the pipeline.

### 1. Pre-requisites:
- Docker Desktop
- NVIDIA Drivers (for CNN Re-scoring and faster MM-GBSA calculations)

### 2. Setup:
Navigate to the *`/docker_variant`* directory and build the Docker image.

**Bash:**
```bash
cd docker_variant
sudo docker build --platform linux/amd64 -t vistadock_p .
```

The pipeline has a GUI that can be used to run it. The GUI Launcher is packaged within the Docker image, so if required, the `gui_launcher.py` script can be extracted from the image using the following command.

**Bash:**
```bash
sudo docker run --rm vistadock_p cat /app/gui_launcher.py > gui_launcher.py
```

### 3. Running the Pipeline:
Open the GUI Launcher using the following command.
```python
python3 gui_launcher.py
```
Read the Disclaimers within the GUI, input all your values in the provided fields, then click on the `GENERATE COMMAND` button. Once the command has been generated in the text box below the button, click on the `RUN PIPELINE` button, and the pipeline will automatically run on your terminal. 




