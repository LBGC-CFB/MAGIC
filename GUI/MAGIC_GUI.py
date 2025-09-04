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
import time
import threading

# REST constants and functions
HEADERS = {"Content-Type": "application/json"}
CUSTOM_CERT = False

# Global optimized cache
_cache = {}
_cache_timestamps = {}
CACHE_DURATION = 3600  # 1 hour

def cached_api_call(cache_key: str, api_function, *args, **kwargs):
    current_time = time.time()
    if (cache_key in _cache and 
        cache_key in _cache_timestamps and 
        current_time - _cache_timestamps[cache_key] < CACHE_DURATION):
        return _cache[cache_key]
    try:
        result = api_function(*args, **kwargs)
        _cache[cache_key] = result
        _cache_timestamps[cache_key] = current_time
        return result
    except Exception as e:
        print(f"Cache error for {cache_key}: {e}")
        return None

def safe_api_call(func, *args, max_retries=2, **kwargs):
    for attempt in range(max_retries):
        try:
            return func(*args, **kwargs)
        except requests.RequestException as e:
            if attempt == max_retries - 1:
                raise e
            time.sleep(1.5 ** attempt)
        except Exception as e:
            if attempt == max_retries - 1:
                raise e
            time.sleep(0.5)
    return None

def get_ensembl_gene_id(gene_name, base_url):
    url = f"{base_url}/xrefs/symbol/homo_sapiens/{gene_name}"
    resp = requests.get(url, headers=HEADERS, verify=CUSTOM_CERT, timeout=10)
    if not resp.ok:
        raise ValueError(f"Error retrieving gene {gene_name}: {resp.text}")
    for e in resp.json():
        if e.get("type") == "gene":
            return e["id"]
    raise ValueError(f"No Ensembl gene ID found for {gene_name}")

def get_gene_complete_info(gene_name, base_url):
    cache_key = f"complete_gene_{gene_name}_{base_url}"

    def fetch_complete_info():
        gene_id = get_ensembl_gene_id(gene_name, base_url)
        url = f"{base_url}/lookup/id/{gene_id}"
        resp = requests.get(url, headers=HEADERS, params={"expand": 1}, verify=CUSTOM_CERT, timeout=10)
        if not resp.ok:
            raise ValueError(f"Error retrieving gene data for {gene_name}: {resp.text}")
        data = resp.json()
        data['gene_id'] = gene_id
        return data
    
    return cached_api_call(cache_key, fetch_complete_info)

def get_mane_transcripts_for_gene(gene_id, base_url):
    cache_key = f"mane_{gene_id}_{base_url}"
    
    def fetch_mane_transcripts():
        try:
            url = f"{base_url}/overlap/id/{gene_id}?feature=mane"
            resp = requests.get(url, headers=HEADERS, verify=CUSTOM_CERT, timeout=10)
            if resp.ok:
                data = resp.json()
                return [item.get('id') for item in data if 'id' in item]
        except Exception as e:
            print(f"MANE API unavailable for {gene_id}: {e}")
        return []
    
    return cached_api_call(cache_key, fetch_mane_transcripts)

def get_xrefs_for_transcript(transcript_id, base_url):
    cache_key = f"xrefs_{transcript_id}_{base_url}"
    
    def fetch_xrefs():
        url = f"{base_url}/xrefs/id/{transcript_id}"
        resp = requests.get(url, headers=HEADERS, verify=CUSTOM_CERT, timeout=10)
        return resp.json() if resp.ok else []
    
    return cached_api_call(cache_key, fetch_xrefs)

def parse_refseq_version(refseq_id):
    if '.' in refseq_id:
        base_id, version = refseq_id.rsplit('.', 1)
        return base_id, version
    return refseq_id, None

def parse_ensembl_transcript_version(enst_id):
    if '.' in enst_id and enst_id.startswith("ENST"):
        base_id, version = enst_id.rsplit('.', 1)
        return base_id, version
    return enst_id, None

def check_refseq_match(xref, refseq_base, requested_version=None):
    if xref.get("dbname") != "RefSeq_mRNA":
        return False, None, None
    
    primary_id = xref.get("primary_id", "")
    display_id = xref.get("display_id", "")
    version = xref.get("version")
    
    if not (primary_id.startswith(refseq_base) or (display_id and display_id.startswith(refseq_base))):
        return False, None, None
    
    if primary_id == refseq_base and version is not None:
        version_str = str(version)
        return True, f"{primary_id}.{version_str}", version_str
    
    if display_id and '.' in display_id:
        base, vers_str = display_id.rsplit('.', 1)
        if base == refseq_base:
            return True, display_id, vers_str
    
    if primary_id == refseq_base:
        return True, primary_id, None
    
    return False, None, None

def is_mane_transcript(transcript, mane_transcript_ids):
    transcript_id = transcript.get('id')
    if not transcript_id:
        return False
    
    if mane_transcript_ids and transcript_id in mane_transcript_ids:
        return True
    
    return transcript.get('gencode_primary') == 1

def get_transcript_exons(transcript_id, base_url):
    cache_key = f"exons_{transcript_id}_{base_url}"
    
    def fetch_exons():
        url = f"{base_url}/lookup/id/{transcript_id}"
        resp = requests.get(url, headers=HEADERS, params={"expand":1}, verify=CUSTOM_CERT, timeout=10)
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
    
    return cached_api_call(cache_key, fetch_exons)

def get_exon_sequence(species, chrom, start, end, strand, base_url):
    cache_key = f"seq_{species}_{chrom}_{start}_{end}_{strand}_{base_url}"
    
    def fetch_sequence():
        region = f"{chrom}:{start}..{end}:{strand}"
        url = f"{base_url}/sequence/region/{species}/{region}"
        resp = requests.get(url, headers={"Content-Type":"text/plain"}, verify=CUSTOM_CERT, timeout=10)
        if not resp.ok:
            raise ValueError(f"Error retrieving sequence: {resp.text}")
        return resp.text
    
    return cached_api_call(cache_key, fetch_sequence)

def validate_transcript_for_gene(transcript_id, gene_name, base_url):
    original_transcript_id = transcript_id.strip()
    print(f"🔍 Validating transcript {original_transcript_id} for gene {gene_name}")
    
    try:
        gene_data = get_gene_complete_info(gene_name, base_url)
        if not gene_data or 'Transcript' not in gene_data:
            raise ValueError(f"Gene {gene_name} not found or has no transcripts")
        
        gene_id = gene_data.get('id') or gene_data.get('gene_id')
        transcripts = gene_data['Transcript']
        mane_transcript_ids = get_mane_transcripts_for_gene(gene_id, base_url)
        
        # Case 1: RefSeq ID (NM_XXXXX.X)
        if original_transcript_id.startswith("NM_"):
            print(f"📝 Processing RefSeq ID: {original_transcript_id}")
            
            base_refseq, requested_version = parse_refseq_version(original_transcript_id)
            
            # Search for ANY matching RefSeq transcript in this gene
            candidates = []
            exact_version_match = None
            mane_candidate = None
            
            for transcript in transcripts:
                enstid = transcript.get("id")
                if not enstid:
                    continue
                
                xrefs = get_xrefs_for_transcript(enstid, base_url)
                for xref in xrefs:
                    is_match, full_refseq, version = check_refseq_match(xref, base_refseq)
                    
                    if is_match:
                        is_mane = is_mane_transcript(transcript, mane_transcript_ids)
                        candidate_info = {
                            'ensembl_id': enstid,
                            'refseq_full': full_refseq,
                            'refseq_base': base_refseq,
                            'refseq_version': version,
                            'ensembl_version': transcript.get('version', 1),
                            'is_mane': is_mane,
                            'requested_id': original_transcript_id
                        }
                        
                        candidates.append(candidate_info)
                        
                        # Check for exact version match
                        if requested_version and version == requested_version:
                            exact_version_match = candidate_info
                        
                        # Check for MANE
                        if is_mane and not mane_candidate:
                            mane_candidate = candidate_info
            
            if not candidates:
                raise ValueError(f"Transcript {original_transcript_id} not found in gene {gene_name}")
            
            # Select best candidate: exact version > MANE > latest version
            if exact_version_match:
                selected_transcript = exact_version_match
                version_found = True
                print(f"✅ Exact version match found: {selected_transcript['refseq_full']}")
            elif mane_candidate:
                selected_transcript = mane_candidate
                version_found = False
                print(f"⚠️ Requested version {requested_version} not found, using MANE transcript: {selected_transcript['refseq_full']}")
            else:
                # Use latest version (highest version number)
                candidates.sort(key=lambda x: int(x['refseq_version']) if x['refseq_version'] and x['refseq_version'].isdigit() else 0, reverse=True)
                selected_transcript = candidates[0]
                version_found = False
                print(f"⚠️ Requested version {requested_version} not found, using latest version: {selected_transcript['refseq_full']}")
            
            return selected_transcript['ensembl_id'], selected_transcript, version_found
        
        # Case 2: Ensembl ID (ENST00000XXXXX.X)
        elif original_transcript_id.startswith("ENST"):
            print(f"📝 Processing Ensembl ID: {original_transcript_id}")
            
            base_tid, requested_version = parse_ensembl_transcript_version(original_transcript_id)
            
            # FIRST: Check if transcript base ID belongs to this gene
            valid_transcript_ids = [t.get('id') for t in transcripts if t.get('id')]
            
            if base_tid not in valid_transcript_ids:
                valid_list = ", ".join(valid_transcript_ids[:5])
                if len(valid_transcript_ids) > 5:
                    valid_list += f" ... (+{len(valid_transcript_ids) - 5} others)"
                raise ValueError(f"Transcript {original_transcript_id} does not belong to gene {gene_name}.\nValid transcripts for {gene_name}: {valid_list}")
            
            print(f"✅ Transcript {base_tid} belongs to gene {gene_name}")
            
            # Get transcript details
            target_transcript = next((t for t in transcripts if t.get('id') == base_tid), None)
            if not target_transcript:
                raise ValueError(f"Transcript {original_transcript_id} not found in gene {gene_name}")
            
            actual_version = target_transcript.get('version', 1)
            is_mane = is_mane_transcript(target_transcript, mane_transcript_ids)
            
            # Check version match and apply fallback if needed
            version_found = True
            if requested_version and str(actual_version) != str(requested_version):
                # Version mismatch - check if we need fallback
                print(f"⚠️ Version {requested_version} requested, but version {actual_version} found for {base_tid}")
                
                # Try to find MANE transcript or use current one
                if not is_mane:
                    # Look for MANE alternative in this gene
                    mane_alternatives = [t for t in transcripts if is_mane_transcript(t, mane_transcript_ids)]
                    if mane_alternatives:
                        target_transcript = mane_alternatives[0]
                        actual_version = target_transcript.get('version', 1)
                        base_tid = target_transcript.get('id')
                        is_mane = True
                        print(f"⚠️ Switched to MANE transcript: {base_tid}.{actual_version}")
                
                version_found = False
            else:
                print(f"✅ Ensembl transcript validated: {base_tid}.{actual_version}")
            
            # Create transcript_info
            transcript_info = {
                'ensembl_id': base_tid,
                'refseq_full': f"{base_tid}.{actual_version}",
                'refseq_base': base_tid,
                'refseq_version': str(actual_version),
                'ensembl_version': actual_version,
                'is_mane': is_mane,
                'version_found': version_found,
                'requested_id': original_transcript_id
            }
            
            return base_tid, transcript_info, version_found
        
        else:
            raise ValueError(f"Unsupported transcript ID format: {original_transcript_id}")
    
    except Exception as e:
        print(f"❌ Validation failed for {original_transcript_id}: {e}")
        raise e

def fetch_exon_sequence(transcript_id, exon_number, gene_name, base_url, species="human"):
    print(f"🔍 Fetching exon {exon_number} for {transcript_id} (gene: {gene_name})")
    
    try:
        # Unified validation for both transcript types
        ensembl_transcript_id, transcript_info, version_found = validate_transcript_for_gene(
            transcript_id, gene_name, base_url
        )
        
        # Retrieve exons
        print(f"📥 Fetching exons for {ensembl_transcript_id}")
        exons, strand = get_transcript_exons(ensembl_transcript_id, base_url)
        
        if not exons:
            raise ValueError(f"No exons found for transcript {ensembl_transcript_id}")
        print(f"📊 Found {len(exons)} exons")
        
        if exon_number < 1 or exon_number > len(exons):
            raise ValueError(f"Exon #{exon_number} out of range (1–{len(exons)}) for {gene_name} ({ensembl_transcript_id})")
        
        # Retrieve sequence
        ex = exons[exon_number-1]
        print(f"🧬 Fetching sequence for exon {exon_number}: {ex['seq_region_name']}:{ex['start']}-{ex['end']}")
        
        sequence = get_exon_sequence(species, ex["seq_region_name"], ex["start"], ex["end"], strand, base_url)
        
        if not sequence:
            raise ValueError(f"Empty sequence returned for exon {exon_number}")
        
        print(f"✅ Sequence retrieved: {len(sequence)} bp")
        return sequence, transcript_info
        
    except Exception as e:
        print(f"❌ Error in fetch_exon_sequence: {e}")
        raise e

class MAGIC_GUI(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.protocol("WM_DELETE_WINDOW", self.on_window_close)

        # Define path of ico and ppm
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
        icon_path = os.path.join(base_path, "MAGIC.ico")
        ppm_path = os.path.join(base_path, "MAGIC.ppm")

        # Define main window
        self.title("MAGIC, create reference files")
        self.iconbitmap(icon_path)
        self.geometry('1500x900')
        self.resizable(True, True)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        # Define canvas and frame (used for scrollbar)
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

        # Define logo and attach to previously created frame
        self.logo = PhotoImage(file=ppm_path)
        self.logo_label = Label(self.scrollable_frame, image=self.logo)
        self.logo_label.grid(row=0, column=1, columnspan=6, pady=20)

        # Define tkinter variables
        self.outdir = StringVar()
        self.seq = StringVar()
        self.UTR5p = StringVar()
        self.UTR3p = StringVar()
        self.gene_name = StringVar()
        self.transcript_name = StringVar()
        self.construction_name = StringVar()
        self.genome_version = StringVar(value="GRCh38")
        self.exon_number = IntVar()
        self.exs_list = []
        self.used_transcript_info = None
        self.loading_widgets = {}

        # Define fonts
        self.ask_font = tkinter.font.Font(family="Helvetica", size="10", weight='bold')
        self.answer_font = tkinter.font.Font(family="Helvetica", size="10")
        self.entry_font = tkinter.font.Font(family="Helvetica", size="10")
        self.button_font = tkinter.font.Font(family="Helvetica", size="10")
        self.button_exit_font = tkinter.font.Font(family="Helvetica", size="17", weight='bold')

        # Build GUI
        self.add_widgets()
        self.values_not_completed = []
        
    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

    def on_window_close(self):
        self.destroy()

    def show_loading(self, widget_key, message="Loading...", row=None):
        if widget_key in self.loading_widgets:
            return
        if hasattr(self, 'exon_inputs_frame'):
            loading_frame = ttk.Frame(self.exon_inputs_frame)
            if row is not None:
                loading_frame.grid(row=row, column=9, padx=5, pady=2, sticky="w")
            else:
                loading_frame.grid(row=0, column=9, padx=5, pady=2, sticky="w")
            progress = ttk.Progressbar(loading_frame, mode='indeterminate', length=100)
            progress.pack(side=tk.LEFT, padx=(0, 5))
            progress.start(interval=50)
            label = ttk.Label(loading_frame, text=message, font=("Helvetica", 9))
            label.pack(side=tk.LEFT)
            self.loading_widgets[widget_key] = loading_frame

    def hide_loading(self, widget_key):
        if widget_key in self.loading_widgets:
            self.loading_widgets[widget_key].destroy()
            del self.loading_widgets[widget_key]

    def show_error(self, error_message):
        messagebox.showerror("Error", error_message)

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
        try:
            gen = generate_files(self)
            if gen.values_not_completed or not gen.success:
                return
            if sys.platform == 'win32':
                os.startfile(os.path.normpath(self.outdir.get()))
            elif sys.platform == 'darwin':
                os.system(f'open "{os.path.normpath(self.outdir.get())}"')
            else:
                os.system(f'xdg-open "{os.path.normpath(self.outdir.get())}"')
        except Exception as e:
            return

    def add_widgets(self):
        # Widget 1: ask output directory
        self.ask_dir = Label(self.scrollable_frame, text='Select output directory', justify="left", font=self.ask_font)
        self.ask_dir.grid(row=5, column=0, padx=5, pady=10, ipadx=3, ipady=3, sticky="w")
        self.ask_btn_dir = Button(self.scrollable_frame, text='Choose directory', command=self.choose_outdir, font=self.button_font)
        self.ask_btn_dir.grid(row=5, column=1, ipadx=1, ipady=1)

        # Widget 2: ask fasta file
        self.ask_input = Label(self.scrollable_frame, text='Select input file', justify="left", font=self.ask_font)
        self.ask_input.grid(row=6, column=0, padx=5, pady=10, ipadx=3, ipady=3, sticky="w")
        self.ask_btn_input = Button(self.scrollable_frame, text='Choose input file', command=self.choose_fasta_file, font=self.button_font)
        self.ask_btn_input.grid(row=6, column=1, ipadx=1, ipady=1)

        # Widget 3: ask exon-specific vector sequences 5' and 3'
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

        # Widget 4, 5 & 6: ask gene, transcript and construction name
        self.ask_gene = Label(self.scrollable_frame, text=' Enter: Gene symbol', justify="left", font=self.ask_font)
        self.ask_gene.grid(row=8, column=0, padx=5, pady=10, sticky="w")
        self.entry_gene = Entry(self.scrollable_frame, font=self.answer_font)
        self.entry_gene.grid(row=8, column=1, padx=5, pady=10, ipadx=3, ipady=3)

        self.ask_tr = Label(self.scrollable_frame, text='Transcript ID\n(NM_XXXXX.X or ENST00000XXXXX.X)', justify="left", font=self.ask_font)
        self.ask_tr.grid(row=8, column=2, padx=5, pady=10, ipadx=3, ipady=3)
        self.entry_tr = Entry(self.scrollable_frame, font=self.answer_font)
        self.entry_tr.grid(row=8, column=3, padx=5, pady=10, ipadx=3, ipady=3)

        self.ask_const = Label(self.scrollable_frame, text='Construction name', justify="left", font=self.ask_font)
        self.ask_const.grid(row=8, column=4, padx=5, pady=10, ipadx=3, ipady=3)
        self.entry_const = Entry(self.scrollable_frame, font=self.answer_font)
        self.entry_const.grid(row=8, column=5, padx=5, pady=10, ipadx=3, ipady=3)

        self.ask_btn_gene_tr_const = Button(self.scrollable_frame, text='Validate ✅', command=self.get_construction_informations, bg="green", fg="white", font=self.ask_font)
        self.ask_btn_gene_tr_const.grid(row=8, column=6, ipadx=1, ipady=1, sticky="e")

        # Widget 7: ask human genome version 
        self.ask_genome_version = Label(self.scrollable_frame, text='Select human genome version', font=self.ask_font)
        self.ask_genome_version.grid(row=9, column=2, padx=5, pady=10, sticky="w")

        self.radio_grch38 = tk.Radiobutton(self.scrollable_frame, text="GRCh38 (default)", variable=self.genome_version, value="GRCh38", font=self.answer_font)
        self.radio_grch38.grid(row=9, column=3, padx=5, pady=5, sticky="w")

        self.radio_grch37 = tk.Radiobutton(self.scrollable_frame, text="GRCh37", variable=self.genome_version, value="GRCh37", font=self.answer_font)
        self.radio_grch37.grid(row=9, column=4, padx=5, pady=5, sticky="w")

        # Widget 8: ask exon(s) number(s)
        self.ask_ex = Label(self.scrollable_frame, text='Select number of exons', font=self.ask_font)
        self.ask_ex.grid(row=9, column=0, padx=5, pady=10, ipadx=3, ipady=3)
        self.ask_btn_ex = Button(self.scrollable_frame, text='Select number of exons', command=self.ask_exon, font=self.button_font)
        self.ask_btn_ex.grid(row=9, column=1, ipadx=1, ipady=1)

        # Final exit button to generate files
        self.exit_button = Button(self.scrollable_frame, text="★⋆⭒˚.⋆  Generate files  ⋆⭒˚.⋆★", command=self.close_open_folder, font=self.button_exit_font, bg="#EC0C6E", fg="#F8D7E5")
        self.exit_button.grid(row=1000, column=2, pady=20, columnspan=3)

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
                with open(self.file_path, "r") as filin:
                    content = filin.read()
                    self.seq.set(content.replace("\n", ""))
            tk.Label(self.scrollable_frame, text=f"sequence length: {len(self.seq.get())} bp", font=self.answer_font, fg="#8D064B").grid(row=6, column=5, padx=5, pady=10)
            self.ask_btn_input.config(relief=SUNKEN, fg="gray")

    def get_UTR_informations(self):
        utr5p = self.entry_seq_ex_5p.get("1.0", "end").strip()
        utr3p = self.entry_seq_ex_3p.get("1.0", "end").strip()
        if not utr5p or not utr3p:
            messagebox.showwarning("Warning", "Please enter both 5' and 3' sequences for constitutive exon vector.")
            return
        self.UTR5p.set(utr5p)
        self.UTR3p.set(utr3p)
        self.ask_btn_seq_ex_5p3p.config(relief=SUNKEN,bg="white", fg="green")

    def get_construction_informations(self):
        gene = self.entry_gene.get().strip()
        tr = self.entry_tr.get().strip()
        const = self.entry_const.get().strip()
        if not gene or not tr or not const:
            messagebox.showwarning("Warning", "Please enter gene symbol, Transcript ID, and Construction name.")
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
            eentry = Entry(self.exon_inputs_frame, font=self.answer_font)
            eentry.grid(row=row, column=2,  padx=5, pady=5, ipadx=3, ipady=3)
            btn = Button(self.exon_inputs_frame, text="Fulfill 🔍", font=self.ask_font, bg="#0077cc", fg="white",
                         command=lambda idx=i: self.fill_single_exon_async(idx))
            btn.grid(row=row, column=3, padx=5, pady=5)
            st = ScrolledText(self.exon_inputs_frame, font=self.answer_font, height=2, width=60)
            st.grid(row=row, column=4,  padx=5, pady=5, ipadx=3, ipady=3, columnspan=4)
            if i < len(previous):
                eentry.insert(0, previous[i][0])
                st.insert("1.0", previous[i][1])
            self.ex_infos.append([eentry, btn, st])
        
        self.ask_ex_seq_batch = Button(self.exon_inputs_frame, text="Fulfill all 🚀", command=self.fill_multiple_exons_batch, bg="#0077cc", fg="white", font=self.ask_font)
        self.ask_ex_seq_batch.grid(row=nb, column=9, padx=5, pady=10, sticky="w")
        
        self.ask_ex_seq = Button(self.exon_inputs_frame, text="Validate ✅", command=lambda:self.keep_exons_info(lst=self.exs_list), bg="green", fg="white", font=self.ask_font)
        self.ask_ex_seq.grid(row=nb, column=10, padx=5, pady=10, sticky="w")
        self.ask_btn_ex.config(relief="raised", state="normal")
        self.activate_right_click_on_all_inputs()

    def fill_single_exon_async(self, idx):
        def background_task():            
            try:
                entry, btn, txt = self.ex_infos[idx]
                exon_num = entry.get().strip()
                
                if not exon_num or not exon_num.isdigit():
                    raise ValueError("Enter a valid exon number")
                
                gene = self.entry_gene.get().strip()
                tr = self.entry_tr.get().strip()
                base_url = self.get_ensembl_base_url()
                
                if not gene or not tr:
                    raise ValueError("Gene name and transcript ID are required")
                
                seq, transcript_info = fetch_exon_sequence(
                    transcript_id=tr, 
                    exon_number=int(exon_num), 
                    gene_name=gene, 
                    base_url=base_url, 
                    species="human"
                )
                
                self.after(0, lambda s=seq, t=transcript_info, i=idx: self.update_exon_result(i, s, t))
                
            except Exception as e:
                error_msg = str(e)
                self.after(0, lambda msg=error_msg: self.show_error(msg))
            finally:
                self.after(0, lambda i=idx: self.hide_loading(f"exon_{i}"))
                self.after(0, lambda i=idx: self.restore_button(i))
        
        entry, btn, txt = self.ex_infos[idx]
        btn.config(state="disabled", text="⏳ Loading...", bg="gray")
        self.show_loading(f"exon_{idx}", f"Exon {idx+1}...", row=idx)
        
        thread = threading.Thread(target=background_task, daemon=True)
        thread.start()

    def restore_button(self, idx):
        if idx < len(self.ex_infos):
            entry, btn, txt = self.ex_infos[idx]
            btn.config(state="normal", text="Fulfill 🔍", bg="#0077cc")

    def update_batch_progress(self, message):
        if "batch" in self.loading_widgets:
            loading_frame = self.loading_widgets["batch"]
            for child in loading_frame.winfo_children():
                if isinstance(child, ttk.Label):
                    child.config(text=message)
                    break
    
    def fill_multiple_exons_batch(self):
        def batch_process():
            gene = self.entry_gene.get().strip()
            tr = self.entry_tr.get().strip()
            base_url = self.get_ensembl_base_url()
            
            try:
                ensembl_transcript_id, transcript_info, version_found = validate_transcript_for_gene(
                    tr, gene, base_url
                )
                
                exons, strand = get_transcript_exons(ensembl_transcript_id, base_url)
                
                processed_count = 0
                total_exons = sum(1 for entry, btn, txt in self.ex_infos 
                                if entry.get().strip() and entry.get().strip().isdigit())
                
                for i, (entry, btn, txt) in enumerate(self.ex_infos):
                    exon_num = entry.get().strip()
                    if exon_num and exon_num.isdigit():
                        try:
                            ex = exons[int(exon_num)-1]
                            seq = get_exon_sequence("human", ex["seq_region_name"], 
                                                  ex["start"], ex["end"], strand, base_url)
                            self.after(0, lambda i=i, s=seq: self.update_single_exon(i, s))
                            processed_count += 1
                            progress_msg = f"Exon {processed_count}/{total_exons}..."
                            self.after(0, self.update_batch_progress, progress_msg)
                        except Exception as e:
                            print(f"Error with exon {exon_num}: {e}")
                
                self.after(0, lambda: setattr(self, 'used_transcript_info', transcript_info))
                self.after(0, self.display_transcript_info)
                
            except Exception as e:
                self.after(0, self.show_error, str(e))
            finally:
                self.after(0, self.hide_loading, "batch")
        
        self.show_loading("batch", "Retrieving all exons sequences...", row=0)
        thread = threading.Thread(target=batch_process, daemon=True)
        thread.start()

    def update_single_exon(self, idx, seq):
        if idx < len(self.ex_infos):
            entry, btn, txt = self.ex_infos[idx]
            txt.delete("1.0", "end")
            txt.insert("1.0", seq)

    def update_exon_result(self, idx, seq, transcript_info):
        entry, btn, txt = self.ex_infos[idx]
        if transcript_info:
            self.used_transcript_info = transcript_info
        txt.delete("1.0", "end")
        txt.insert("1.0", seq)
        self.display_transcript_info()

    def display_transcript_info(self):
        if not self.used_transcript_info:
            return
            
        info = self.used_transcript_info
        
        if hasattr(self, 'transcript_info_label'):
            self.transcript_info_label.destroy()
        
        requested_transcript = self.entry_tr.get().strip()
        is_ensembl_input = requested_transcript.startswith("ENST")
        
        if is_ensembl_input:
            info_text = f"Used transcript:"
            info_text += f"\n\tEnsembl: {info['ensembl_id']}"
            if info.get('ensembl_version'):
                info_text += f".{info['ensembl_version']}"
            
            if 'version_found' in info and not info['version_found']:
                base_requested, version_requested = parse_ensembl_transcript_version(requested_transcript)
                info_text += f"\n⚠️  Transcript version {version_requested} not found \nin Ensembl REST API"
                info_text += f"\nSwitched to {info['refseq_full']}"
        else:
            info_text = f"Used transcript:"
            info_text += f"\n\tRefSeq: {info['refseq_full']}"
            info_text += f"\n\tEnsembl: {info['ensembl_id']}"
            if info.get('ensembl_version'):
                info_text += f".{info['ensembl_version']}"
            
            if '.' in requested_transcript and info['refseq_full'] != requested_transcript:
                info_text += f"\n⚠️  Transcript version {requested_transcript} not found \nin Ensembl REST API"
                info_text += f"\nSwitched to {info['refseq_full']}"
        
        if info.get('is_mane'):
            info_text += "\n\tType: MANE Select ⭐"
        
        self.transcript_info_label = Label(
            self.exon_inputs_frame, 
            text=info_text, 
            font=self.answer_font, 
            fg="#0077cc", 
            justify="left",
            relief="solid",
            borderwidth=1,
            padx=5,
            pady=5,
            bg="#f0f8ff"
        )
        
        max_row = len(self.ex_infos) if hasattr(self, 'ex_infos') else 0
        self.transcript_info_label.grid(row=0, column=11, rowspan=max_row+1, padx=10, pady=5, sticky="n")

    def keep_exons_info(self, lst):
        lst.clear()
        errors = []
        for i, ex_inf in enumerate(self.ex_infos):
            tmp_lst = []
            exon_number = ex_inf[0].get().strip()
            exon_seq = ex_inf[2].get("1.0", "end-1c").strip()
            if not exon_number:
                errors.append(f"Exon {i+1}: missing exon number.")
            if not exon_seq:
                errors.append(f"Exon {i+1}: missing exon sequence.")
            tmp_lst.append(exon_number)
            tmp_lst.append(exon_seq)
            lst.append(tmp_lst)

        if errors:
            messagebox.showerror("Validation Error", "\n".join(errors))
        else:
            tk.Label(self.exon_inputs_frame, text=f"{len(self.ex_infos)} exon(s) and sequence(s) validated", font=self.answer_font, fg="#8D064B").grid(row=len(self.ex_infos), column=5, padx=5, pady=10, sticky="w")
            self.ask_ex_seq.config(relief=SUNKEN, bg="white", fg="green")

class generate_files():
    def __init__(self, gui):
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

        self.values = {"output directory": self.outdir, "input file": self.seq, "constitutive exon vector 5 prime sequence": self.UTR5p, "constitutive exon vector 3 prime sequence": self.UTR3p, "gene symbol": self.gene_name, "transcript name": self.transcript_name, "construction name": self.construction_name, "exon number and sequence": self.exs_list}
        self.values_not_completed = [val_name for val_name, value in self.values.items() if not value]

        if self.values_not_completed:
            message = "These fields are not defined: \n\n\t" + "\n\t".join(self.values_not_completed) + "\n\nPlease verify and click 'Validate'"
            messagebox.showerror("ERROR. Files not generated", message)
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
            filout.write(f'#MAGIC creation of artificial genome fasta file for {self.construction_name} construction used for minigene splicing assays and massively parallel sequencing\n')
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
                filout.write(f'#MAGIC creation of artificial genome gtf file for {self.construction_name} construction used for minigene splicing assays and massively parallel sequencing\n')
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
