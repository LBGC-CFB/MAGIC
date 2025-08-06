import os
import sys
import tkinter as tk
import tkinter.font
from tkinter import ttk
from tkinter import messagebox
from tkinter import Label, Entry, Button, StringVar, IntVar, PhotoImage, SUNKEN
from tkinter.filedialog import askdirectory, askopenfilename
from tkinter.scrolledtext import ScrolledText
import requests



# REST constants and functions
HEADERS = {"Content-Type": "application/json"}
CUSTOM_CERT = False

def get_ensembl_gene_id(gene_name, base_url):
    url = f"{base_url}/xrefs/symbol/homo_sapiens/{gene_name}"
    resp = requests.get(url, headers=HEADERS, verify=CUSTOM_CERT)
    if not resp.ok:
        raise ValueError(f"Error retrieving gene {gene_name}: {resp.text}")
    for e in resp.json():
        if e.get("type") == "gene":
            return e["id"]
    raise ValueError(f"No Ensembl gene ID found for {gene_name}")

def get_transcripts_for_gene(gene_id, base_url):
    url = f"{base_url}/lookup/id/{gene_id}"
    resp = requests.get(url, headers=HEADERS, params={"expand":1}, verify=CUSTOM_CERT)
    if not resp.ok:
        raise ValueError(f"Error retrieving transcripts for {gene_id}: {resp.text}")
    return resp.json().get("Transcript", [])

def get_xrefs_for_transcript(transcript_id, base_url):
    url = f"{base_url}/xrefs/id/{transcript_id}"
    resp = requests.get(url, headers=HEADERS, verify=CUSTOM_CERT)
    if not resp.ok:
        return []
    return resp.json()

def find_ensembl_transcript(gene_name, refseq_id, base_url):
    try:
        gene_id = get_ensembl_gene_id(gene_name, base_url)
        trans = get_transcripts_for_gene(gene_id, base_url)
        for t in trans:
            enstid = t.get("id")
            if not enstid:
                continue
            x = get_xrefs_for_transcript(enstid, base_url)
            if any(xref["dbname"]=="RefSeq_mRNA" and xref["primary_id"]==refseq_id for xref in x):
                return enstid
    except ValueError:
        return None
    return None

def get_transcript_exons(transcript_id, base_url):
    url = f"{base_url}/lookup/id/{transcript_id}"
    resp = requests.get(url, headers=HEADERS, params={"expand":1}, verify=CUSTOM_CERT)
    if not resp.ok:
        raise ValueError(f"Error retrieving transcript {transcript_id}: {resp.text}")
    data = resp.json()
    ex = data.get("Exon", [])
    strand = data.get("strand")
    if ex and "exon_number" in ex[0]:
        sorted_ex = sorted(ex, key=lambda x: x["exon_number"], reverse=(strand==-1))
    else:
        sorted_ex = sorted(ex, key=lambda x: x["start"], reverse=(strand==-1))
    return sorted_ex, strand

def get_exon_sequence(species, chrom, start, end, strand, base_url):
    region = f"{chrom}:{start}..{end}:{strand}"
    url = f"{base_url}/sequence/region/{species}/{region}"
    print(f"url: {url}")
    resp = requests.get(url, headers={"Content-Type":"text/plain"}, verify=CUSTOM_CERT)
    if not resp.ok:
        raise ValueError(f"Error retrieving sequence: {resp.text}")
    return resp.text

def fetch_exon_sequence(transcript_id, exon_number, gene_name, base_url, species="human"):
    tid = transcript_id.strip()
    if tid.startswith("NM_"):
        mt = find_ensembl_transcript(gene_name, tid, base_url)
        if mt:
            tid = mt
        else:
            raise ValueError(f"Could not find {tid} for gene {gene_name}. Please verify NM.")
    exons, strand = get_transcript_exons(tid, base_url)
    if exon_number < 1 or exon_number > len(exons):
        raise ValueError(f"Exon #{exon_number} out of range (1–{len(exons)}) for {gene_name} ({tid})")
    ex = exons[exon_number-1]
    return get_exon_sequence(species, ex["seq_region_name"], ex["start"], ex["end"], strand, base_url)


class MAGIC_GUI(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.protocol("WM_DELETE_WINDOW", self.on_window_close)

        # define path of ico and ppm
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
        #icon_path = os.path.join(base_path, "MAGIC.xbm") #linux
        icon_path = os.path.join(base_path, "MAGIC.ico")
        ppm_path = os.path.join(base_path, "MAGIC.ppm")

        # define main window
        self.title("MAGIC, create reference files")
        #self.iconbitmap("@"+icon_path) #linux
        self.iconbitmap(icon_path)
        self.geometry('1500x900')
        self.resizable(True, True)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        # define canvas and frame (used for scrollbar)
        self.frame = ttk.Frame(self)
        self.frame.grid(row=0, column=0, sticky="nsew")

        self.canvas = tk.Canvas(self.frame)
        self.scrollbarv = ttk.Scrollbar(self.frame, orient="vertical", command=self.canvas.yview)
        self.scrollbarh = ttk.Scrollbar(self.frame, orient="horizontal", command=self.canvas.xview)

        self.scrollable_frame = ttk.Frame(self.canvas)
        self.scrollable_frame.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))

        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=self.scrollbarv.set)
        self.canvas.configure(xscrollcommand=self.scrollbarh.set)

        self.scrollbarv.grid(row=0, column=1, sticky="ns")
        self.scrollbarh.grid(row=1, column=0, sticky="ew")
        self.canvas.grid(row=0, column=0, sticky="nsew")

        self.frame.grid(row=0, column=0, sticky="nsew")
        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(0, weight=1)
        self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)

        # define logo and attach to previously created frame
        self.logo = PhotoImage(file=ppm_path)
        self.logo_label = Label(self.scrollable_frame, image=self.logo)
        self.logo_label.grid(row=0, column=1, columnspan=6, pady=20)

        # define tkinter variables
        self.outdir = StringVar()
        self.seq = StringVar()
        self.UTR5p = StringVar()
        self.UTR3p = StringVar()
        self.gene_name = StringVar()
        self.transcript_name = StringVar()
        self.construction_name = StringVar()
        self.genome_version = StringVar(value="GRCh38")  # Default value
        self.exon_number = IntVar()
        self.exs_list = []

        # define fonts
        self.ask_font = tkinter.font.Font(family="Helvetica", size="10", weight='bold')
        self.answer_font = tkinter.font.Font(family="Helvetica", size="10")
        self.entry_font = tkinter.font.Font(family="Helvetica", size="10")
        self.button_font = tkinter.font.Font(family="Helvetica", size="10")
        self.button_exit_font = tkinter.font.Font(family="Helvetica", size="17", weight='bold')

        # build gui
        self.add_widgets()

        # To store missing fields after validation
        self.values_not_completed = []

    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

    def on_window_close(self):
        self.destroy()

    def add_right_click_menu(self, widget):
        menu = tk.Menu(widget, tearoff=0)
        menu.add_command(label="Cut", command=lambda: widget.event_generate("<<Cut>>"))
        menu.add_command(label="Copy", command=lambda: widget.event_generate("<<Copy>>"))
        menu.add_command(label="Paste", command=lambda: widget.event_generate("<<Paste>>"))

        def show_menu(event):
            try:
                menu.tk_popup(event.x_root, event.y_root)
            finally:
                menu.grab_release()

        widget.bind("<Button-3>", show_menu)
        widget.bind("<Control-Button-1>", show_menu)
    
    def activate_right_click_on_all_inputs(self):
        def recursive_attach(widget):
            for child in widget.winfo_children():
                if isinstance(child, (tk.Entry, ScrolledText)):
                    self.add_right_click_menu(child)
                recursive_attach(child)
        recursive_attach(self)

    def close_open_folder(self):
        # validate fields and generate files
        try:
            gen = generate_files(self)
            if gen.values_not_completed or not gen.success:
                # error message shown in generate_files
                return
            # If files created successfully, open folder and close app
            if sys.platform == 'win32':
                os.startfile(os.path.normpath(self.outdir.get()))
            elif sys.platform == 'darwin':
                os.system(f'open "{os.path.normpath(self.outdir.get())}"')
            else:
                os.system(f'xdg-open "{os.path.normpath(self.outdir.get())}"')
            self.destroy()
        except Exception as e:
            return

    def add_widgets(self):
        # widget1: ask output directory
        self.ask_dir = Label(self.scrollable_frame, text='Select ouput directory', justify="left", font=self.ask_font)
        self.ask_dir.grid(row=5, column=0, padx=5, pady=10, ipadx=3, ipady=3, sticky="w")
        self.ask_btn_dir = Button(self.scrollable_frame, text='Choose directory', command=self.choose_outdir, font=self.button_font)
        self.ask_btn_dir.grid(row=5, column=1, ipadx=1, ipady=1)

        # widget2: ask fasta file
        self.ask_input = Label(self.scrollable_frame, text='Select input file', justify="left", font=self.ask_font)
        self.ask_input.grid(row=6, column=0, padx=5, pady=10, ipadx=3, ipady=3, sticky="w")
        self.ask_btn_input = Button(self.scrollable_frame, text='Choose input file', command=self.choose_fasta_file, font=self.button_font)
        self.ask_btn_input.grid(row=6, column=1, ipadx=1, ipady=1)

        # widget3: ask exon-specific vector sequences 5' and 3'
        self.ask_seq_ex_5p = Label(self.scrollable_frame, text='Enter: \nconstitutive exon vector \n5 prime sequence', justify="left", font=self.ask_font)
        self.ask_seq_ex_5p.grid(row=7, column=0, padx=5, pady=10, sticky="w")
        self.entry_seq_ex_5p = ScrolledText(self.scrollable_frame, font=self.answer_font, height=1, width=35)
        self.entry_seq_ex_5p.grid(row=7, column=1, padx=5, pady=10, ipadx=3, ipady=3)

        self.ask_seq_ex_3p = Label(self.scrollable_frame, text='constitutive exon vector \n3 prime sequence', justify="left", font=self.ask_font)
        self.ask_seq_ex_3p.grid(row=7, column=2, padx=5, pady=10, sticky="w")
        self.entry_seq_ex_3p = ScrolledText(self.scrollable_frame, font=self.answer_font, height=1, width=35)
        self.entry_seq_ex_3p.grid(row=7, column=3, padx=5, pady=10, ipadx=3, ipady=3)

        self.ask_btn_seq_ex_5p3p = Button(self.scrollable_frame, text='Validate ✅', command=self.get_UTR_informations, bg="green", fg="white", font=self.ask_font)
        self.ask_btn_seq_ex_5p3p.grid(row=7, column=4, ipadx=1, ipady=1, sticky="e")

        # widget4, 5 & 6: ask gene, transcript and construction name
        self.ask_gene = Label(self.scrollable_frame, text=' Enter: Gene name', justify="left", font=self.ask_font)
        self.ask_gene.grid(row=8, column=0, padx=5, pady=10, sticky="w")
        self.entry_gene = Entry(self.scrollable_frame, font=self.answer_font)
        self.entry_gene.grid(row=8, column=1, padx=5, pady=10, ipadx=3, ipady=3)

        self.ask_tr = Label(self.scrollable_frame, text='Transcript name', justify="left", font=self.ask_font)
        self.ask_tr.grid(row=8, column=2, padx=5, pady=10, ipadx=3, ipady=3)
        self.entry_tr = Entry(self.scrollable_frame, font=self.answer_font)
        self.entry_tr.grid(row=8, column=3, padx=5, pady=10, ipadx=3, ipady=3)

        self.ask_const = Label(self.scrollable_frame, text='Construction name', justify="left", font=self.ask_font)
        self.ask_const.grid(row=8, column=4, padx=5, pady=10, ipadx=3, ipady=3)
        self.entry_const = Entry(self.scrollable_frame, font=self.answer_font)
        self.entry_const.grid(row=8, column=5, padx=5, pady=10, ipadx=3, ipady=3)

        self.ask_btn_gene_tr_const = Button(self.scrollable_frame, text='Validate ✅', command=self.get_construction_informations, bg="green", fg="white", font=self.ask_font)
        self.ask_btn_gene_tr_const.grid(row=8, column=6, ipadx=1, ipady=1, sticky="e")

        # widget 7: ask human genome version 
        self.ask_genome_version = Label(self.scrollable_frame, text='Select human genome version', font=self.ask_font)
        self.ask_genome_version.grid(row=9, column=2, padx=5, pady=10, sticky="w")

        self.radio_grch38 = tk.Radiobutton(self.scrollable_frame, text="GRCh38 (default)", variable=self.genome_version, value="GRCh38", font=self.answer_font)
        self.radio_grch38.grid(row=9, column=3, padx=5, pady=5, sticky="w")

        self.radio_grch37 = tk.Radiobutton(self.scrollable_frame, text="GRCh37", variable=self.genome_version, value="GRCh37", font=self.answer_font)
        self.radio_grch37.grid(row=9, column=4, padx=5, pady=5, sticky="w")

        # widget 8: ask exon(s) number(s)
        self.ask_ex = Label(self.scrollable_frame, text='Select number of exons', font=self.ask_font)
        self.ask_ex.grid(row=9, column=0, padx=5, pady=10, ipadx=3, ipady=3)
        self.ask_btn_ex = Button(self.scrollable_frame, text='Select number of exons', command=self.ask_exon, font=self.button_font)
        self.ask_btn_ex.grid(row=9, column=1, ipadx=1, ipady=1)

        # final exit button to generate files
        self.exit_button = Button(self.scrollable_frame, text="★⋆⭒˚.⋆  Generate files  ⋆⭒˚.⋆★", command=self.close_open_folder, font=self.button_exit_font, bg="#EC0C6E", fg="#F8D7E5")
        self.exit_button.grid(row=1000, column=2, pady=20, columnspan=3)

        #enables mouse clicks
        self.activate_right_click_on_all_inputs()

    def choose_outdir(self):
        self.selected_dir = askdirectory()
        if self.selected_dir:
            tk.Label(self.scrollable_frame, text=f"{self.selected_dir}", font=self.answer_font, fg="#8D064B").grid(row=5, column=2, padx=5, pady=10, columnspan=4, sticky="w")
            self.outdir.set(self.selected_dir)
            self.ask_btn_dir.config(relief=SUNKEN, fg="gray")

    def choose_fasta_file(self):
        self.file_path = askopenfilename(filetypes=[("Fasta files", "*.fa *.fasta"), ("GenBank files", "*.gb"), ("All files", "*.*")])
        if self.file_path:
            tk.Label(self.scrollable_frame, text=f"{self.file_path}", font=self.answer_font, fg="#8D064B").grid(row=6, column=2, padx=5, pady=10, columnspan=4, sticky="w")
            ext = self.file_path.split(".")[-1].lower()
            seq_list = []
            if ext == "gb":
                with open(self.file_path, "r") as filin:
                    for line in filin:
                        line = line.strip().split()
                        if line and line[0].isdigit():
                            seq_list.extend(line[1:])
                self.seq.set("".join(seq_list))
            elif ext in ("fasta", "fa"):
                with open(self.file_path, "r") as filin:
                    seq = []
                    for line in filin:
                        if not line.startswith(">"):
                            seq.append(line.strip())
                self.seq.set("".join(seq))
            else:
                # fallback for unknown file types - try to read entire file as sequence
                with open(self.file_path, "r") as filin:
                    content = filin.read()
                    self.seq.set(content.replace("\n", ""))
            tk.Label(self.scrollable_frame, text=f"sequence length: {len(self.seq.get())} bp", font=self.answer_font, fg="#8D064B").grid(row=6, column=5, padx=5, pady=10)
            self.ask_btn_input.config(relief=SUNKEN, fg="gray")

    def get_UTR_informations(self):
        # Store sequences 5p and 3p for constitutive exon vector
        utr5p = self.entry_seq_ex_5p.get("1.0", "end").strip()
        utr3p = self.entry_seq_ex_3p.get("1.0", "end").strip()
        if not utr5p or not utr3p:
            messagebox.showwarning("Warning", "Please enter both 5' and 3' sequences for constitutive exon vector.")
            return
        self.UTR5p.set(utr5p)
        self.UTR3p.set(utr3p)
        self.ask_btn_seq_ex_5p3p.config(relief=SUNKEN,bg="white", fg="green")

    def get_construction_informations(self):
        # Store gene, transcript and construction names
        gene = self.entry_gene.get().strip()
        tr = self.entry_tr.get().strip()
        const = self.entry_const.get().strip()
        if not gene or not tr or not const:
            messagebox.showwarning("Warning", "Please enter gene name, transcript name, and construction name.")
            return
        self.gene_name.set(gene)
        self.transcript_name.set(tr)
        self.construction_name.set(const)
        self.ask_btn_gene_tr_const.config(relief=SUNKEN, bg="white", fg="green")

    def get_ensembl_base_url(self):
        if self.genome_version.get() == "GRCh37":
            return "https://grch37.rest.ensembl.org"
        else:
            return "https://rest.ensembl.org"

    def ask_exon(self):
        nb = tk.simpledialog.askinteger("Number of exons?", "Number of exons?", parent=self)
        if nb is None:
            return
        self.exon_number.set(nb)
        previous = []
        if hasattr(self, 'ex_infos'):
            for e, btn, txt in self.ex_infos:
                previous.append((e.get().strip(), txt.get("1.0","end-1c").strip()))
        if hasattr(self, 'exon_inputs_frame'):
            self.exon_inputs_frame.destroy()
        self.exon_inputs_frame = ttk.Frame(self.scrollable_frame)
        self.exon_inputs_frame.grid(row=11, column=0, columnspan=10, padx=10, pady=10, sticky="w")
        self.ex_infos = []
        for i in range(nb):
            row = i
            Label(self.exon_inputs_frame, text=f'Exon {i+1} number?', font=self.ask_font).grid(row=row, column=1, padx=5, pady=5, ipadx=3, ipady=3, sticky="w")
            eentry = Entry(self.exon_inputs_frame, font=self.answer_font)#, width=10)
            eentry.grid(row=row, column=2,  padx=5, pady=5, ipadx=3, ipady=3)#padx=5, pady=5)
            btn = Button(self.exon_inputs_frame, text="Fullfill 🔍", font=self.ask_font, bg="#0077cc", fg="white",
                         command=lambda idx=i: self.fill_single_exon(idx))
            btn.grid(row=row, column=3, padx=5, pady=5)
            st = ScrolledText(self.exon_inputs_frame, font=self.answer_font, height=2, width=60)
            st.grid(row=row, column=4,  padx=5, pady=5, ipadx=3, ipady=3, columnspan=4)
            if i < len(previous):
                eentry.insert(0, previous[i][0])
                st.insert("1.0", previous[i][1])
            self.ex_infos.append([eentry, btn, st])
        self.ask_ex_seq = Button(self.exon_inputs_frame, text="Validate ✅", command=lambda:self.keep_exons_info(lst=self.exs_list), bg="green", fg="white", font=self.ask_font)
        self.ask_ex_seq.grid(row=nb, column=10, padx=5, pady=10, sticky="w")
        self.ask_btn_ex.config(relief="raised", state="normal")
        self.activate_right_click_on_all_inputs()

    def fill_single_exon(self, idx):
        try:
            entry, btn, txt = self.ex_infos[idx]
            exon_num = entry.get().strip()
            if not exon_num or not exon_num.isdigit():
                raise ValueError("Enter a valid exon number")
            gene = self.entry_gene.get().strip()
            tr = self.entry_tr.get().strip()
            base_url = self.get_ensembl_base_url()
            seq = fetch_exon_sequence(tr, int(exon_num), gene, base_url=base_url)
            txt.delete("1.0", "end")
            txt.insert("1.0", seq)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def keep_exons_info(self, lst):
        lst.clear()
        erreurs = []
        for i, ex_inf in enumerate(self.ex_infos):
            tmp_lst = []
            exon_number = ex_inf[0].get().strip()
            exon_seq = ex_inf[2].get("1.0", "end-1c").strip()
            if not exon_number:
                erreurs.append(f"Exon {i+1}: missing exon number.")
            if not exon_seq:
                erreurs.append(f"Exon {i+1}: missing exon sequence.")
            tmp_lst.append(exon_number)
            tmp_lst.append(exon_seq)
            lst.append(tmp_lst)

        if erreurs:
            messagebox.showerror("Validation Error", "\n".join(erreurs))
        else:
            tk.Label(self.exon_inputs_frame, text=f"{len(self.ex_infos)} exon(s) and sequence(s) validated", font=self.answer_font, fg="#8D064B").grid(row=len(self.ex_infos), column=5, padx=5, pady=10, sticky="w")
            self.ask_ex_seq.config(relief=SUNKEN, bg="white", fg="green")



class generate_files():
    def __init__(self, gui):
        # retrieve value from MAGIC GUI
        self.gui = gui
        self.gene_name = magic.gene_name.get().strip()
        self.outdir = magic.outdir.get().strip()
        self.seq = magic.seq.get().strip()
        self.UTR5p = magic.UTR5p.get().strip()
        self.UTR3p = magic.UTR3p.get().strip()
        self.transcript_name = magic.transcript_name.get().strip()
        self.construction_name = magic.construction_name.get().strip()
        self.exon_number = magic.exon_number.get()
        self.exs_list = magic.exs_list
        self.success = True

        self.values = {"output directory": self.outdir, "input file": self.seq, "constitutive exon vector 5 prime sequence": self.UTR5p, "constitutive exon vector 3 prime sequence": self.UTR3p, "gene name": self.gene_name, "transcript name": self.transcript_name, "construction name": self.construction_name, "exon number and sequence": self.exs_list}
        self.values_not_completed = [val_name for val_name, value in self.values.items() if not value]

        # generate files if all fields are completed
        if self.values_not_completed:
            message = "These fields are not defined: \n\n\t" + "\n\t".join(self.values_not_completed) + "\n\nPlease verify and click 'Validate'"
            messagebox.showerror("ERROR. Files not generated", message)
            print("These fields are not defined: " + ", ".join(self.values_not_completed))
            self.success = False
            return
            
        self.create_output_dir()
        if not self.create_fasta_file():
            self.success = False
            return
        if not self.create_gtf_file():
            self.success = False
            return
        if self.success:
            messagebox.showinfo("Success", "All files have been generated successfully!")

    def validate_fields(self):
        # Check all necessary fields filled
        if not self.gui.outdir.get():
            self.values_not_completed.append("Output directory")
        if not self.gui.seq.get():
            self.values_not_completed.append("Input sequence (fasta or genbank)")
        if not self.gui.UTR5p.get():
            self.values_not_completed.append("Constitutive exon vector 5' sequence")
        if not self.gui.UTR3p.get():
            self.values_not_completed.append("Constitutive exon vector 3' sequence")
        if not self.gui.gene_name.get():
            self.values_not_completed.append("Gene name")
        if not self.gui.transcript_name.get():
            self.values_not_completed.append("Transcript name")
        if not self.gui.construction_name.get():
            self.values_not_completed.append("Construction name")
        if self.gui.exon_number.get() == 0:
            self.values_not_completed.append("Number of exons")
        if not self.gui.exs_list or any(not exon[1].strip() for exon in self.gui.exs_list):
            self.values_not_completed.append("Exon sequences")

    def create_output_dir(self):
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def create_fasta_file(self):
        self.seq_search_UTR5p_start = self.seq.lower().find(self.UTR5p.lower())
        self.seq_search_UTR3p_end = self.seq.lower().find(self.UTR3p.lower())

        if self.seq_search_UTR5p_start == -1:
            messagebox.showwarning("Warning", "5' sequence not found in input file. Please verify sequence and input file.")
            return False
        if self.seq_search_UTR3p_end == -1:
            messagebox.showwarning("Warning", "3' sequence not found in input file. Please verify sequence and input file.")
            return False

        self.insert_start = self.seq[self.seq_search_UTR5p_start+len(self.UTR5p):]
        self.insert_end = self.seq[:self.seq_search_UTR3p_end]
        self.seq_split = self.UTR5p + self.insert_start + self.insert_end + self.UTR3p
        
        with open(self.outdir+f"/{self.construction_name}.fasta", "w") as filout:
            filout.write(f'#MAGIC creation of artficial genome fasta file for {self.construction_name} construction used for minigene splicing assays and massively parallel sequencing\n')
            filout.write(f">{self.construction_name}\n")
            filout.write(self.seq_split)
        return True

    def create_gtf_file(self):
        seq_UTR5_start = self.seq_split.lower().find(self.UTR5p.lower()); seq_UTR5_end = seq_UTR5_start + len(self.UTR5p)
        seq_UTR3_start = self.seq_split.lower().find(self.UTR3p.lower()); seq_UTR3_end = seq_UTR3_start + len(self.UTR3p)
        strand = "+"
        
        if seq_UTR5_start == -1:
            messagebox.showwarning("Warning", "5' sequence not found in input file. Please verify sequence and input file.")
            return False
        if seq_UTR3_start == -1:
            messagebox.showwarning("Warning", "3' sequence not found in input file. Please verify sequence and input file.")
            return False
        try:
            with open (self.outdir+f"/{self.construction_name}.gtf", "w") as filout:
                filout.write(f'#MAGIC creation of artficial genome gtf file for {self.construction_name} construction used for minigene splicing assays and massively parallel sequencing\n')
                filout.write(f'{self.construction_name}\tMAGIC\ttranscript\t1\t{len(self.seq_split)}\t.\t{strand}\t.\tgene_name "{self.gene_name}"; transcript_id "{self.transcript_name}_{self.construction_name}";\n')
                filout.write(f'{self.construction_name}\tMAGIC\tUTR\t{str(seq_UTR5_start+1)}\t{seq_UTR5_end}\t.\t{strand}\t.\tgene_name "{self.gene_name}"; transcript_id "{self.transcript_name}_{self.construction_name}";\n')
                for ex in self.exs_list:
                    ex_nb = ex[0].strip(); ex_seq = ex[1].strip()
                    self.seq_search_start = self.seq_split.lower().find(ex_seq.lower())
                    self.seq_search_end = self.seq_search_start + len(ex_seq)
                    if self.seq_search_start == -1:
                        messagebox.showwarning("Warning", f"Exon {ex_nb} sequence not found in input file. Please verify sequence and input file.")
                        return False

                    filout.write(f'{self.construction_name}\tMAGIC\tCDS\t{str(self.seq_search_start+1)}\t{str(self.seq_search_end)}\t.\t{strand}\t.\tgene_name "{self.gene_name}"; transcript_id "{self.transcript_name}_{self.construction_name}"; exon_number "{ex_nb}";\n')
                    filout.write(f'{self.construction_name}\tMAGIC\texon\t{str(self.seq_search_start+1)}\t{str(self.seq_search_end)}\t.\t{strand}\t.\tgene_name "{self.gene_name}"; transcript_id "{self.transcript_name}_{self.construction_name}"; exon_number "{ex_nb}";\n')
                filout.write(f'{self.construction_name}\tMAGIC\tUTR\t{str(seq_UTR3_start+1)}\t{seq_UTR3_end}\t.\t{strand}\t.\tgene_name "{self.gene_name}"; transcript_id "{self.transcript_name}_{self.construction_name}";\n')
            return True
        except Exception as e:
            if os.path.exists(self.outdir+f"/{self.construction_name}.gtf"):
                try:
                    os.remove(self.outdir+f"/{self.construction_name}.gtf")
                except Exception as remove_error:
                    messagebox.showerror("Cleanup Error", f"An error occurred while removing incomplete GTF file:\n{remove_error}")



if __name__ == '__main__':
    magic = MAGIC_GUI()
    magic.mainloop()