import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import multiprocessing
import subprocess
import threading

max_sys_cpus = multiprocessing.cpu_count()

class PipelineGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("VISTADOCK-P Command Constructor")
        self.root.geometry("800x850")

        style = ttk.Style()
        style.configure("Bold.TLabel", font=("Helvetica", 13, 'bold'))
        style.configure("Section.TLabel", font=("Helvetica", 12, 'bold')) 

        main_frame = tk.Frame(root)
        main_frame.pack(fill=tk.BOTH, expand=1)

        self.canvas = tk.Canvas(main_frame)
        self.scrollbar = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=self.canvas.yview)
        self.scrollable_frame = ttk.Frame(self.canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )

        self.window_id = self.canvas.create_window((0,0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=self.scrollbar.set)
        
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfig(self.window_id, width=e.width))

        self.canvas.pack(side="left", fill="both", expand=True)
        self.scrollbar.pack(side="right", fill="y")

        self.create_widgets()

    def create_widgets(self):
        paddings = {'padx':10, 'pady':5}

        disclaimer_frame = ttk.LabelFrame(self.scrollable_frame, text="IMPORTANT DISCLAIMERS")
        disclaimer_frame.pack(fill="x", padx=10, pady=(10,20))

        warnings = ("1. Ensure that you have built the Docker Image for this tool using the instructions provided on the GitHub Repository.\n"
                    "2. Ensure that Docker is running before starting the pipeline.\n"
                    "3. Ensure that your Receptor PDB is fixed before running the pipeline (bond orders, valencies, non-standard residues etc.).\n"
                    "4. Fields marked with an asterisk (*) are mandatory and the pipeline will not run if they are not filled in.\n"
                    "5. Ensure that your Receptor PDB, Ligand File and Output Directory are in the same directory you are running the pipeline from."
        )

        ttk.Label(
            disclaimer_frame,
            text=warnings,
            style="Section.TLabel",
            justify="left",
        ).pack(anchor="w", padx=10, pady=10)

        ttk.Label(self.scrollable_frame, text="1. File Selection", style="Section.TLabel").pack(anchor="w", pady=(20,5), padx=10)

        self.receptor_path = self.create_file_input("Receptor PDB*: ", "Select Receptor PDB")
        self.ligand_path = self.create_file_input("Ligand(s) (SDF/Mol2)*: ", "Select Ligand File (SDF/Mol2/PDB)")
        self.output_path = self.create_dir_input("Output Directory: ", "pipeline_results")

        ttk.Label(self.scrollable_frame, text="2. Grid Box Definition", style="Section.TLabel").pack(anchor="w", pady=(20,5), padx=10)

        box_frame = ttk.Frame(self.scrollable_frame)
        box_frame.pack(fill='x', padx=10)

        self.cx = self.create_labelled_entry(box_frame, "Center X*", 0, 0)
        self.cy = self.create_labelled_entry(box_frame, "Center Y*", 0, 2)
        self.cz = self.create_labelled_entry(box_frame, "Center Z*", 0, 4)

        self.sx = self.create_labelled_entry(box_frame, "Size X*", 1, 0, default = "22")
        self.sy = self.create_labelled_entry(box_frame, "Size Y*", 1, 2, default = "22")
        self.sz = self.create_labelled_entry(box_frame, "Size Z*", 1, 4, default = "22")

        ttk.Label(self.scrollable_frame, text="3. Docking Parameters", style="Section.TLabel").pack(anchor="w", pady=(20,5), padx=10)

        self.dock_frame = ttk.Frame(self.scrollable_frame)
        self.dock_frame.pack(fill='x', padx=10)

        self.exh_rapid = self.create_labelled_entry(self.dock_frame, "Exhaustiveness (RAPID): ", 0, 0, default="2")
        self.exh_balanced = self.create_labelled_entry(self.dock_frame, "Exhaustiveness (BALANCED): ", 1, 0, default="8")
        self.exh_ultra = self.create_labelled_entry(self.dock_frame, "Exhaustiveness (ULTRA): ", 2, 0, default="32")

        self.frac_rapid = self.create_labelled_entry(self.dock_frame, "Fraction to keep after RAPID: ", 0, 2, default="0.5")
        self.frac_balanced = self.create_labelled_entry(self.dock_frame, "Fraction to keep after BALANCED: ", 1, 2, default="0.3")
        self.frac_ultra = self.create_labelled_entry(self.dock_frame, "Fraction to keep after ULTRA: ", 2, 2, default="1.0")

        self.single_step_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(self.dock_frame, text="SingleStep Docking (Bypass RAPID & BALANCED, run only ULTRA)", variable=self.single_step_var, command=self.toggle_single_step).grid(row=3, column=0, columnspan=4, sticky='w', pady=(10, 0))

        ttk.Label(self.scrollable_frame, text="4. Workflow Options", style="Section.TLabel").pack(anchor="w", pady=(20,5), padx=10)

        opt_frame = ttk.Frame(self.scrollable_frame)
        opt_frame.pack(fill='x', padx=10)

        ttk.Label(opt_frame, text=f"CPU Limit (Docking) (0=Max) ({max_sys_cpus} available): ").grid(row=0, column=0, sticky='w', pady=5)
        self.cpu_var = tk.StringVar(value="0")
        self.cpu_entry = tk.Entry(opt_frame, textvariable=self.cpu_var, width=10)
        self.cpu_entry.grid(row=0, column=1, sticky='w', padx=5)
        self.cpu_var.trace_add("write", lambda *args: self.validate_cpu_input(self.cpu_var, self.cpu_entry))

        ttk.Label(opt_frame, text=f"CPU Limit (Ligand Prep) (0=Auto) ({max_sys_cpus} available): ").grid(row=1, column=0, sticky='w', pady=5)
        self.prep_cpu_var = tk.StringVar(value="0")
        self.prep_cpu_entry = tk.Entry(opt_frame, textvariable=self.prep_cpu_var, width=10)
        self.prep_cpu_entry.grid(row=1, column=1, sticky='w', padx=5)
        self.prep_cpu_var.trace_add("write", lambda *args: self.validate_cpu_input(self.prep_cpu_var, self.prep_cpu_entry))

        ttk.Label(opt_frame, text="Force Field: ").grid(row=2, column=0, sticky='w', pady=5)
        self.ff_var = tk.StringVar(value="amber14-all")
        ff_opts = ["amber14-all",
                   "amber14-SB",
                   "amber19-SB",
                   "amber99sb",
                   "amber99sbildn",
                   "amber03",
                   "amber10",
                   "amber96",
                   "charmm36",
                   "amoeba2013",
                   "amoeba2018"
                   ]
        ttk.OptionMenu(opt_frame, self.ff_var, "amber14-all", *ff_opts).grid(row=2, column=1, sticky='w', padx=5)


        ttk.Label(opt_frame, text="Scoring Mode: ").grid(row=3, column=0, sticky='w', pady=5)
        self.scoring_var = tk.StringVar(value="both")
        scoring_opts = ["both", "vina", "vinardo"]
        ttk.OptionMenu(opt_frame, self.scoring_var, "both", *scoring_opts).grid(row=3, column=1, sticky='w', padx=5)

        self.lipinksi_var = tk.BooleanVar(value=False)
        self.no_mmgbsa_var = tk.BooleanVar(value=False)
        self.no_plip_var = tk.BooleanVar(value=False)
        
        self.run_prep_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(opt_frame, text="Perform Ligand Preparation (Uncheck to bypass)", variable=self.run_prep_var, command=self.toggle_prep_options).grid(row=4, column=0, columnspan=2, sticky='w', pady=(10,2))
        
        self.prep_opts_frame = ttk.Frame(opt_frame)
        self.prep_opts_frame.grid(row=5, column=0, columnspan=2, sticky='w', padx=20)
        
        ttk.Label(self.prep_opts_frame, text="Ligand Preparation Force Field: ").grid(row=0, column=0, sticky='w')
        self.prep_ff_var = tk.StringVar(value="MMFF94")
        prep_ffs = ["MMFF94", "MMFF94s", "UFF", "Ghemical", "GAFF"]
        ttk.OptionMenu(self.prep_opts_frame, self.prep_ff_var, "MMFF94", *prep_ffs).grid(row=0, column=1, sticky='w', padx=5)

        ttk.Checkbutton(opt_frame, text="Apply Lipinski's Rule of 5 Filter", variable=self.lipinksi_var).grid(row=6, column=0, columnspan=2, sticky='w')
        ttk.Checkbutton(opt_frame, text="Skip MM-GBSA (Faster but less accurate)", variable=self.no_mmgbsa_var).grid(row=7, column=0, columnspan=2, sticky='w')
        ttk.Checkbutton(opt_frame, text="Skip PLIP Analysis (Faster but no interaction data)", variable=self.no_plip_var).grid(row=8, column=0, columnspan=2, sticky='w')

        self.cnn_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(opt_frame, text="Enable GNINA CNN Re-Scoring (Requires GPU)", variable=self.cnn_var).grid(row=9, column=0, columnspan=2, sticky='w')

        ttk.Label(opt_frame, text="GPU Device ID (0,1, ...): ").grid(row=10, column=0, sticky='w', pady=5)
        self.gpu_id_var = tk.StringVar(value="0")
        tk.Entry(opt_frame, textvariable=self.gpu_id_var, width=5).grid(row=10, column=1, sticky='w', padx=5)

        ttk.Label(self.scrollable_frame, text="5. MM-GBSA Parameters", style="Section.TLabel").pack(anchor='w', pady=(20,5), padx=10)
        md_frame = ttk.Frame(self.scrollable_frame)
        md_frame.pack(fill='x', padx=10)

        self.temp_var = self.create_labelled_entry(md_frame, "Temperature (K): ", 0, 0, default="300.0")
        self.fric_var = self.create_labelled_entry(md_frame, "Co-efficient of Friction (1/ps): ", 2, 0, default="1.0")
        self.step_var = self.create_labelled_entry(md_frame, "Timestep (femtoseconds): ", 4, 0, default="2.0")

        btn_frame = ttk.Frame(self.scrollable_frame)
        btn_frame.pack(pady=20)

        gen_btn = ttk.Button(btn_frame, text="GENERATE COMMAND", command=self.generate_command)
        gen_btn.pack(side='left', padx=10, ipadx=10, ipady=5)

        run_btn = ttk.Button(btn_frame, text="RUN PIPELINE", command=self.run_pipeline)
        run_btn.pack(side='left', padx=10, ipadx=10, ipady=5)

        ttk.Label(self.scrollable_frame, text="Copy this command into your terminal: ", style="Bold.TLabel").pack(anchor='w', padx=10)
        self.output_text = tk.Text(self.scrollable_frame, height=12, width=90)
        self.output_text.pack(padx=10, pady=5, fill='both', expand=True)

    def create_file_input(self, label_text, placeholder):
        frame = ttk.Frame(self.scrollable_frame)
        frame.pack(fill='x', padx=10, pady=2)
        ttk.Label(frame, text=label_text, width=25, anchor='w').pack(side='left')

        entry_var = tk.StringVar()
        entry = ttk.Entry(frame, textvariable=entry_var)

        def browse():
            filename = filedialog.askopenfilename()
            if filename:  entry_var.set(filename)

        ttk.Button(frame, text="Browse...", command=browse).pack(side='right')
        entry.pack(side='left', fill='x', expand=True, padx=5)
        return entry_var
    
    def create_dir_input(self, label_text, default):
        frame = ttk.Frame(self.scrollable_frame)
        frame.pack(fill='x', padx=10, pady=2)
        ttk.Label(frame, text=label_text, width=25, anchor='w').pack(side='left')

        entry_var = tk.StringVar()
        entry = ttk.Entry(frame, textvariable=entry_var)

        def browse():
            dirname = filedialog.askdirectory()
            if dirname: entry_var.set(dirname)

        ttk.Button(frame, text="Browse...", command=browse).pack(side='right')
        entry.pack(side='left', fill='x', expand=True, padx=5)
        return entry_var
    
    def create_labelled_entry(self, parent, text, row, col, default=""):
        ttk.Label(parent, text=text).grid(row=row, column=col, sticky='e', padx=5, pady=2)
        entry_var = tk.StringVar(value=default)
        ttk.Entry(parent, textvariable=entry_var, width=10).grid(row=row, column=col+1, sticky='w', padx=5, pady=2)
        return entry_var
    
    def toggle_prep_options(self):
        state = 'normal' if self.run_prep_var.get() else 'disabled'
        for child in self.prep_opts_frame.winfo_children():
            child.configure(state=state)
            
    def toggle_single_step(self):
        state = 'disabled' if self.single_step_var.get() else 'normal'
        for row in [0,1]:
            for child in self.dock_frame.grid_slaves(row=row):
                child.configure(state=state)
    
    def generate_command(self):
        rec_full_path = self.receptor_path.get()
        rec = rec_full_path.split('/')
        rec = f"/data/{rec[-1]}"
        lig_full_path = self.ligand_path.get()
        lig = lig_full_path.split('/')
        lig = f"/data/{lig[-1]}"
        output_full_path = self.output_path.get()
        out = output_full_path.split('/')
        out = f"/data/{out[-1]}"

        if not rec or not lig:
            messagebox.showerror("Missing Files", "Please select both Receptor and Ligand files.")
            return
        
        cmd_parts = ["sudo docker run --rm --gpus all -v $(pwd):/data -u $(id -u):$(id -g) vistadock_p python3 -u vistadockp.py "]
        
        cmd_parts.append(f'--receptor_pdb "{rec}"')
        cmd_parts.append(f'--ligand "{lig}"')
        cmd_parts.append(f'--output_dir "{out}"')

        try:
            cmd_parts.append(f'--center_x {self.cx.get()} --center_y {self.cy.get()} --center_z {self.cz.get()}')
            cmd_parts.append(f'--size_x {self.sx.get()} --size_y {self.sy.get()} --size_z {self.sz.get()}')
        except:
            messagebox.showerror("Error", "Box coordinates and dimensions must be numbers.")
            return
        
        cmd_parts.append(f'--exh_rapid {self.exh_rapid.get()} --exh_balanced {self.exh_balanced.get()} --exh_ultra {self.exh_ultra.get()}')
        cmd_parts.append(f'--frac_rapid {self.frac_rapid.get()} --frac_balanced {self.frac_balanced.get()} --frac_ultra {self.frac_ultra.get()}')

        cmd_parts.append(f'--scoring {self.scoring_var.get()}')

        try:
            cmd_parts.extend(['--cpu_count', self.cpu_var.get()])
        except:
            cmd_parts.extend(['--cpu_count', '0'])
        
        if not self.run_prep_var.get():
            cmd_parts.append('--skip_prep')
        else:
            cmd_parts.append(f'--prep_ff {self.prep_ff_var.get()}')
        
        try:
            cmd_parts.extend(['--prep_cpus', self.prep_cpu_var.get()])
        except:
            cmd_parts.extend(['--prep_cpus', '0'])

        cmd_parts.append(f'--forcefield {self.ff_var.get()}')

        if self.lipinksi_var.get(): cmd_parts.append('--lipinski')
        if self.single_step_var.get(): cmd_parts.append('--single_step')
        if self.no_mmgbsa_var.get(): cmd_parts.append('--no_mmgbsa')
        if self.no_plip_var.get(): cmd_parts.append('--no_plip')
        if self.cnn_var.get(): cmd_parts.append('--cnn_scoring rescore')
        cmd_parts.append(f'--gpu_device {self.gpu_id_var.get()}')

        cmd_parts.append(f'--gbsa_temp {self.temp_var.get()}')
        cmd_parts.append(f'--gbsa_friction {self.fric_var.get()}')
        cmd_parts.append(f'--gbsa_timestep {self.step_var.get()}')

        full_command = " ".join(cmd_parts)
        self.output_text.delete(1.0, tk.END)
        self.output_text.insert(tk.END, full_command)

    def run_pipeline(self):
        """Runs the command generated by the user directly in the terminal."""
        cmd_str = self.output_text.get("1.0", tk.END).strip()

        if not cmd_str:
            messagebox.showwarning("Empty Command", "Please generate a command first!")
            return
        
        if not messagebox.askyesno("Confirm Run?", "This will start the pipeline in your terminal.\nThe GUI will remain active.\n\nProceed?"):
            return
        
        t = threading.Thread(target=self._execute_subprocess, args=(cmd_str,))
        t.start()

    def _execute_subprocess(self, cmd_str):
        """Worker thread that runs the command."""
        try:
            print("\n" + "="*60)
            print(">>> STARTING PIPELINE")
            print("="*60 + "\n")
            subprocess.run(cmd_str, check=True, shell=True)

            messagebox.showinfo("Success", "Pipeline run finished successfully!")

        except subprocess.CalledProcessError as e:
            messagebox.showerror("Pipeline Failed", f"The pipeline encountered an error.\nExit Code: {e.returncode}\nCheck Terminal for details.")
        except FileNotFoundError:
            messagebox.showerror("Error", "Could not find 'python3'. Are you in the correct enviornment?")
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred:\n{e}")

    def validate_cpu_input(self, var, entry_widget):
        """Checks if input is > max_sys_cpus and turns background red."""
        value = var.get()

        entry_widget.config(bg="white")

        if value.isdigit():
            if int(value) > max_sys_cpus:
                entry_widget.config(bg="#ce6161")
        elif value == "":
            pass
        else:
            entry_widget.config(bg="#ce6161")

if __name__ == "__main__":
    root = tk.Tk()

    root.lift()
    root.attributes('-topmost', True)
    root.after_idle(root.attributes, '-topmost', False)

    app = PipelineGUI(root)
    root.mainloop()