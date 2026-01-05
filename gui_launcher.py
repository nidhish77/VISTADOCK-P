import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os

class PipelineGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Virtual Screening Pipeline Command Constructor")
        self.root.geometry("700x850")

        style = ttk.Style()
        style.configure("Bold.TLabel", font=("Helvetica", 13, 'bold'))
        style.configure("Section.TLabel", font=("Helvetica", 12, 'bold'), foreground="#333") 

        main_frame = tk.Frame(root)
        main_frame.pack(fill=tk.BOTH, expand=1)

        self.canvas = tk.Canvas(main_frame)
        self.scrollbar = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=self.canvas.yview)
        self.scrollable_frame = ttk.Frame(self.canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )

        self.canvas.create_window((0,0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        self.canvas.pack(side="left", fill="both", expand=True)
        self.scrollbar.pack(side="right", fill="y")

        self.create_widgets()

    def create_widgets(self):
        paddings = {'padx':10, 'pady':5}
        
        ttk.Label(self.scrollable_frame, text="1. File Selection", style="Section.TLabel").pack(anchor="w", pady=(10,5), padx=10)

        self.receptor_path = self.create_file_input("Receptor PDB: ", "Select Receptor PDB")
        self.ligand_path = self.create_file_input("Ligand File(SDF/Mol2): ", "Select Ligand File (SDF/Mol2/PDB)")
        self.output_path = self.create_dir_input("Output Directory: ", "pipeline_results")

        ttk.Label(self.scrollable_frame, text="2. Grid Box Definition", style="Section.TLabel").pack(anchor="w", pady=(10,5), padx=10)

        box_frame = ttk.Frame(self.scrollable_frame)
        box_frame.pack(fill='x', padx=10)

        self.cx = self.create_labelled_entry(box_frame, "Center X", 0, 0)
        self.cy = self.create_labelled_entry(box_frame, "Center Y", 0, 2)
        self.cz = self.create_labelled_entry(box_frame, "Center Z", 0, 4)

        self.sx = self.create_labelled_entry(box_frame, "Size X", 1, 0, default = "22")
        self.sy = self.create_labelled_entry(box_frame, "Size Y", 1, 2, default = "22")
        self.sz = self.create_labelled_entry(box_frame, "Size Z", 1, 4, default = "22")

        ttk.Label(self.scrollable_frame, text="3. Docking Parameters", style="Section.TLabel").pack(anchor="w", pady=(15,5), padx=10)

        dock_frame = ttk.Frame(self.scrollable_frame)
        dock_frame.pack(fill='x', padx=10)

        self.exh_rapid = self.create_labelled_entry(dock_frame, "Exhaustiveness (RAPID): ", 0, 0, default="2")
        self.exh_standard = self.create_labelled_entry(dock_frame, "Exhaustiveness (STANDARD): ", 1, 0, default="8")
        self.exh_precise = self.create_labelled_entry(dock_frame, "Exhaustiveness (PRECISE): ", 2, 0, default="32")

        self.frac_rapid = self.create_labelled_entry(dock_frame, "Fraction to keep after RAPID: ", 0, 2, default="0.5")
        self.frac_standard = self.create_labelled_entry(dock_frame, "Fraction to keep after STANDARD: ", 1, 2, default="0.3")
        self.frac_precise = self.create_labelled_entry(dock_frame, "Fraction to keep after PRECISE: ", 2, 2, default="1.0")

        ttk.Label(self.scrollable_frame, text="4. Workflow Options", style="Section.TLabel").pack(anchor="w", pady=(15,5), padx=10)

        opt_frame = ttk.Frame(self.scrollable_frame)
        opt_frame.pack(fill='x', padx=10)

        ttk.Label(opt_frame, text="Total CPU Limit (Docking) (0=Max): ").grid(row=0, column=0, sticky='w', pady=5)
        self.cpu_var = tk.StringVar(value="0")
        ttk.Entry(opt_frame, textvariable=self.cpu_var, width=10).grid(row=0, column=1, sticky='w', padx=5)

        ttk.Label(opt_frame, text="CPU Limit (Ligand Prep) (0=Auto): ").grid(row=1, column=0, sticky='w', pady=5)
        self.prep_cpu_var = tk.StringVar(value="0")
        ttk.Entry(opt_frame, textvariable=self.prep_cpu_var, width=10).grid(row=1, column=1, sticky='w', padx=5)

        ttk.Label(opt_frame, text="Scoring Mode: ").grid(row=2, column=0, sticky='w', pady=5)
        self.scoring_var = tk.StringVar(value="both")
        scoring_opts = ["both", "vina", "vinardo"]
        ttk.OptionMenu(opt_frame, self.scoring_var, "both", *scoring_opts).grid(row=2, column=1, sticky='w', padx=5)

        self.lipinksi_var = tk.BooleanVar(value=False)
        self.no_mmgbsa_var = tk.BooleanVar(value=False)
        self.no_plip_var = tk.BooleanVar(value=False)

        ttk.Checkbutton(opt_frame, text="Apply Lipinski's Rule of 5 Filter", variable=self.lipinksi_var).grid(row=3, column=0, columnspan=2, sticky='w')
        ttk.Checkbutton(opt_frame, text="Skip MM-GBSA (Faster but less accurate)", variable=self.no_mmgbsa_var).grid(row=4, column=0, columnspan=2, sticky='w')
        ttk.Checkbutton(opt_frame, text="Skip PLIP Analysis (Faster but no interaction data)", variable=self.no_plip_var).grid(row=5, column=0, columnspan=2, sticky='w')

        gen_btn = ttk.Button(self.scrollable_frame,text="GENERATE COMMAND",command=self.generate_command)
        gen_btn.pack(pady=20, ipadx=10, ipady=5)

        ttk.Label(self.scrollable_frame, text="Copy this command into your terminal: ", style="Bold.TLabel").pack(anchor='w', padx=10)
        self.output_text = tk.Text(self.scrollable_frame, height=12, width=90)
        self.output_text.pack(padx=10, pady=5, fill='both', expand=True)

    def create_file_input(self, label_text, placeholder):
        frame = ttk.Frame(self.scrollable_frame)
        frame.pack(fill='x', padx=10, pady=2)
        ttk.Label(frame, text=label_text, width=15, anchor='w').pack(side='left')

        entry_var = tk.StringVar()
        entry = ttk.Entry(frame, textvariable=entry_var)
        entry.pack(side='left', fill='x', expand=True, padx=5)

        def browse():
            filename = filedialog.askopenfilename()
            if filename:  entry_var.set(filename)

        ttk.Button(frame, text="Browse...", command=browse).pack(side='right')
        return entry_var
    
    def create_dir_input(self, label_text, default):
        frame = ttk.Frame(self.scrollable_frame)
        frame.pack(fill='x', padx=10, pady=2)
        ttk.Label(frame, text=label_text, width=15, anchor='w').pack(side='left')

        entry_var = tk.StringVar()
        entry = ttk.Entry(frame, textvariable=entry_var)
        entry.pack(side='left', fill='x', expand=True, padx=5)

        def browse():
            dirname = filedialog.askdirectory()
            if dirname: entry_var.set(dirname)

        ttk.Button(frame, text="Browse...", command=browse).pack(side='right')
        return entry_var
    
    def create_labelled_entry(self, parent, text, row, col, default=""):
        ttk.Label(parent, text=text).grid(row=row, column=col, sticky='e', padx=5, pady=2)
        entry_var = tk.StringVar(value=default)
        ttk.Entry(parent, textvariable=entry_var, width=10).grid(row=row, column=col+1, sticky='w', padx=5, pady=2)
        return entry_var
    
    def generate_command(self):
        rec = self.receptor_path.get()
        lig = self.ligand_path.get()

        if not rec or not lig:
            messagebox.showerror("Missing Files", "Please select both Receptor and Ligand files.")
            return
        
        cmd_parts = ["python3 virtual_screening_pipeline.py "]
        
        cmd_parts.append(f'--receptor_pdb "{rec}"')
        cmd_parts.append(f'--ligand "{lig}"')
        cmd_parts.append(f'--output_dir "{self.output_path.get()}"')

        try:
            cmd_parts.append(f'--center_x {self.cx.get()} --center_y {self.cy.get()} --center_z {self.cz.get()}')
            cmd_parts.append(f'--size_x {self.sx.get()} --size_y {self.sy.get()} --size_z {self.sz.get()}')
        except:
            messagebox.showerror("Error", "Box coordinates and dimensions must be numbers.")
            return
        
        cmd_parts.append(f'--exh_rapid {self.exh_rapid.get()} --exh_standard {self.exh_standard.get()} --exh_precise {self.exh_precise.get()}')
        cmd_parts.append(f'--frac_rapid {self.frac_rapid.get()} --frac_standard {self.frac_standard.get()} --frac_precise {self.frac_precise.get()}')

        cmd_parts.append(f'--scoring {self.scoring_var.get()}')

        try:
            cmd_parts.extend(['--cpu_count', self.cpu_var.get()])
        except:
            cmd_parts.extend(['--cpu_count', '0'])
        
        try:
            cmd_parts.extend(['--prep_cpus', self.prep_cpu_var.get()])
        except:
            cmd_parts.extend(['--prep_cpus', '0'])

        if self.lipinksi_var.get(): cmd_parts.append('--lipinski')
        if self.no_mmgbsa_var.get(): cmd_parts.append('--no_mmgbsa')
        if self.no_plip_var.get(): cmd_parts.append('--no_plip')

        full_command = " ".join(cmd_parts)
        self.output_text.delete(1.0, tk.END)
        self.output_text.insert(tk.END, full_command)

if __name__ == "__main__":
    root = tk.Tk()

    root.lift()
    root.attributes('-topmost', True)
    root.after_idle(root.attributes, '-topmost', False)

    app = PipelineGUI(root)
    root.mainloop()